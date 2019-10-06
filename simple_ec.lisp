;;;
;;; A very simple framework for Evolutionary Computation,
;;; especially genetic programming (GP) and grammar-based genetic programming related
;;; Simplify the previous one, which in retrospect maybe a bit over-architectured.
;;; Basically a bunch of useful routines for custom EC algorithm, and sample routines
;;;   for overall EC flow, that is intended to be copied and modified for custom use.

(defpackage :simple-evolutionary-computation
  (:use :common-lisp)
  (:nicknames :EC)
  (:export :EC :make-individual :individual-p :individual-fitness :individual-chromosome
	   :fill-all-with :fill-with :bernoulli :bernoulli-1/2 :uniform-pick :random-pick :random-choice-c :random-choice
	   :new-individual :higher-fitness-than :higher-fitness-of :lower-fitness-than :lower-fitness-of :best-of
	   :fitness-proportionate-prepare :fitness-proportionate-selector
	   :tournament-prepare :new-2tournament-selector
	   :*crossover-fitness-changes* :*mutation-fitness-changes* :*f-class-of-fitness* :*fitness-record-output*
	   :reset-crossover-fitness-record :reset-mutation-fitness-record :reset-n-fitness-class :report-fitness-changes
	   :add-fitness-record-crossover :add-fitness-record-mutation
	   :new-fitness-record-ind-crossoveror :new-fitness-record-ind-mutator
	   :generational-elitism-EC :*EC-output*))
(in-package :simple-evolutionary-computation)

;;;
;;; Now the simple EC, can be easily modified as needed
(defstruct individual
  "An individual is a problem specific chromosome and a fitness value.
When not evaluated yet, the fitness is NIL. So this provides simple caching to avoid repeated evaluation of the same chromosome."
  fitness
  chromosome)

;;; simple utility functions
(defmacro fill-all-with (v &body body)
  (let ((nv (gensym))
	(ni (gensym)))
    `(let ((,nv ,v))
       (dotimes (,ni (length ,nv))
	 (setf (aref ,nv ,ni) (progn ,@body)))
       ,nv)))
; (fill-all-with (make-array 10)
;   (print "Hello")
;   (random 10))
;;
(defmacro fill-with (v &body parts)
  ;; fill in the vector v with different things at different sections
  ;; parts is list of (n . body) which specifies n slots to fill, and each using (progn body)
  (let ((nv (gensym))
	(ni (gensym))
	(nj (gensym)))
    `(let ((,nv ,v)
	   (,ni 0))
       ,@(mapcar #'(lambda (clause)
		     `(dotimes (,nj ,(car clause))
			(setf (aref ,nv ,ni) (progn ,@(cdr clause)))
			(incf ,ni)))
		 parts)
       ,nv)))
;;
; (fill-with (make-array 10)
;   (2 (+ 1 2))
;   (6 5)
;   (2 (random 1.0)))
;;
(defun bernoulli (p)
  (<= (random 1.0) p))
(defun bernoulli-1/2 ()
  (zerop (random 2)))

(defun uniform-pick (lst)
  (nth (random (length lst)) lst))
(defun random-pick (probs)
  "probs is a vector of positive numbers that add to 1.
  probs represents the distribution of a finite discrete random variable (0 to length(probs)-1).
  Generate a value according to this distribution."
  (let ((r (random 1.0))
	(s 0))
    (min (dotimes (i (length probs))
	   (incf s (aref probs i))
	   (if (< r s)
	       (return i)))
	 (1- (length probs)))))
(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun normalize-probs! (probs)
    (let ((s 0))
      (dotimes (i (length probs))
	(incf s (aref probs i)))
      (dotimes (i (length probs))
	(setf (aref probs i) (/ (aref probs i) s)))
      probs))
  (defun range (start end &optional (step 1))
    (loop for i from start below end by step collect i)))
(defmacro random-choice-c (&rest choices)
  "Each choice is a list where car is probability (may be unnormalized), and cdr is a list of expression (implicitly in progn) to produce a value.
  The probabilities for each choice should be constants."
  (let ((probs (normalize-probs! (coerce (mapcar #'car choices) 'vector))))
    `(case (random-pick ,probs)
       ,@(mapcar (lambda (i c) `((,i) ,@(cdr c)))
		 (range 0 (length choices))
		 choices))))

(defmacro random-choice (&rest choices)
  "Each choice is a list where car is probability, and cdr is a list of expression (implicitly in progn) to produce a value.
  The probabilities for each choice need not be constants. The probability of the last choice is ignored.
  Suitable for small number of choices."
  (let ((s (gensym))
	(r (gensym)))
    (labels ((clause (lst)
	       (if (null (cdr lst))
		   `(progn ,@(cdar lst))
		   `(progn
		      (incf ,s ,(caar lst))
		      (if (< ,r ,s)
			  (progn ,@(cdar lst))
			  ,(clause (cdr lst)))))))
      `(let ((,s 0)
	     (,r (random 1.0)))
	 ,(clause choices)))))
;;(random-choice (0.1 "a") (0.2 "b") (0.4 "c") (0.3 "d"))

;;
(defun new-individual (chr chr-evaluator)
  (make-individual :fitness (funcall chr-evaluator chr)
		   :chromosome chr))
;;
(defun higher-fitness-than (ind1 ind2)
  (> (individual-fitness ind1)
     (individual-fitness ind2)))
(defun higher-fitness-of (ind1 ind2)
  ;; returns the one with higher fitness among ind1 and ind2
  (if (higher-fitness-than ind1 ind2) ind1 ind2))
(defun lower-fitness-than (ind1 ind2)
  (< (individual-fitness ind1)
     (individual-fitness ind2)))
(defun lower-fitness-of (ind1 ind2)
  ;; returns the one with lower fitness among ind1 and ind2
  (if (lower-fitness-than ind1 ind2) ind1 ind2))

(defun best-of (v &optional (predicate #'>))
  (let ((r (aref v 0)))
    (dotimes (i (1- (length v)))
      (if (funcall predicate (aref v (1+ i)) r)
	  (setf r (aref v (1+ i)))))
    r))
; selection
; fitness proportionate
(defun fitness-proportionate-prepare (pop)
  "Gives the preparer for fitness-proportionate selection.
Use the accumulated normalized array of fitness, and the population as the information for selection.
Fitness values are assumed non-negative (better positive) and higher is better."
  (let* ((v (map 'vector #'individual-fitness pop))
	 (total (reduce #'+ v))
	 (acc 0))
    (dotimes (i (length v))
      (setf acc (+ acc (/ (aref v i) total)))
      (setf (aref v i) acc))
    (cons v pop)))

(defun first-smaller-index (x v)
  "Gives the index of vector v which has values smaller than x.
Gives the last index if none is smaller than x."
  (let ((n (length v)))
    (do ((i 0 (1+ i)))
	((or (>= i n) (< (aref v i) x))
	 (if (>= i n) (1- (length v)) i)))))
(defun fitness-proportionate-selector (p)
  ;; p is a pair of the vector of normalized fitness, and the population itself
  (aref (cdr p) (first-smaller-index (random 1.0) (car p))))

; tournament
; Use this to initialize (once at application start) the state for random before use
; (setf *random-state* (make-random-state t))
(defun tournament-prepare (pop) pop)
(defun tournament-random-select (pop)
  (aref pop (random (length pop))))
(defun new-2tournament-selector (better-individual)
  "Gives a deterministic 2 tournament generator."
  #'(lambda (pop)
      (funcall better-individual
	       (tournament-random-select pop)
	       (tournament-random-select pop))))
;;;;;
;; wrapper over chromosome crossoveror to create individual crossoveror that record something related to fitness
(defvar *n-fitness-class* 10)
(defstruct fitness-classes
  total-count   ;; sum of entries in distribution
  distribution  ;; each is count of the corresponding fitness class
  )
(defun new-fitness-classes ()
  (make-fitness-classes :total-count 0 :distribution (make-array *n-fitness-class* :initial-element 0)))
(defun print-fitness-classes (fc-occ &optional (out *standard-output*))
  ;; in fact the default printer for structure looks nice too
  (format out "[total: ~a, counts: ~a]" (fitness-classes-total-count fc-occ) (fitness-classes-distribution fc-occ)))
(defun one-more-case-to (fc-occ c)
  (incf (fitness-classes-total-count fc-occ))
  (incf (aref (fitness-classes-distribution fc-occ) c)))

(defvar *crossover-fitness-changes* nil) ;; element at (i,j) is the fitness-classes of class i crossed with class j
(defvar *mutation-fitness-changes* nil) ;; element at i is fitness-classes of mutating class i
;; should be a unary function to give a class (0 to *n-fitness-class*-1) for a fitness
;; *f-class-of-fitness* should be overided as appropriate
(defvar *f-class-of-fitness* #'(lambda (fitness) (floor fitness 10)))
(defvar *fitness-record-output* *standard-output*)
(defmacro class-of-fitness (f)
  `(funcall *f-class-of-fitness* ,f))
(defun reset-crossover-fitness-record ()
  (setf *crossover-fitness-changes* (make-array (list *n-fitness-class* *n-fitness-class*)))
  (dotimes (i *n-fitness-class*)
    (dotimes (j *n-fitness-class*)
      (setf (aref *crossover-fitness-changes* i j) (new-fitness-classes)))))
(defun reset-mutation-fitness-record ()
  (setf *mutation-fitness-changes* (make-array *n-fitness-class*))
  (fill-all-with *mutation-fitness-changes*
    (new-fitness-classes)))
(defun reset-n-fitness-class (n-class)
  (setf *n-fitness-class* n-class)
  (reset-crossover-fitness-record)
  (reset-mutation-fitness-record))
(defun report-fitness-changes (n)
  ;; report the current statistics for the nth generation
  (format *fitness-record-output* ";;; ============ Generation ~a~%" n)
  (format *fitness-record-output* ";; Crossover (~a classes)~%" *n-fitness-class*)
  (dotimes (i *n-fitness-class*)
    (dotimes (j *n-fitness-class*)
      (format *fitness-record-output* "[~a]x[~a]: " i j)
      (print-fitness-classes (aref *crossover-fitness-changes* i j) *fitness-record-output*)
      (terpri *fitness-record-output*)))
  ;;
  (format *fitness-record-output* "~%;; Mutation (~a classes)~%" *n-fitness-class*)
  (dotimes (i *n-fitness-class*)
    (format *fitness-record-output* "[~a]: " i)
    (print-fitness-classes (aref *mutation-fitness-changes* i) *fitness-record-output*)
    (terpri *fitness-record-output*))
  (format *fitness-record-output* "~%;;; ============================~%"))

(defun add-fitness-record-crossover (ind-crossoveror)
  #'(lambda (ind1 ind2)
      (let ((f1 (individual-fitness ind1))
	    (f2 (individual-fitness ind2))
	    (ind (funcall ind-crossoveror ind1 ind2)))
	(one-more-case-to (aref *crossover-fitness-changes*
				(class-of-fitness f1)
				(class-of-fitness f2))
			  (class-of-fitness (individual-fitness ind)))
	ind)))
(defun new-fitness-record-ind-crossoveror (chr-crossoveror chr-evaluator)
  (add-fitness-record-crossover
   #'(lambda (ind1 ind2)
       (new-individual (funcall chr-crossoveror
				(individual-chromosome ind1)
				(individual-chromosome ind2))
		       chr-evaluator))))
(defun add-fitness-record-mutation (ind-mutator)
  #'(lambda (ind1)
      (let ((f1 (individual-fitness ind1))
	    (ind (funcall ind-mutator ind1)))
	(one-more-case-to (aref *mutation-fitness-changes*
				(class-of-fitness f1))
			  (class-of-fitness (individual-fitness ind)))
	ind)))
(defun new-fitness-record-ind-mutator (chr-mutator chr-evaluator)
  (add-fitness-record-mutation
   #'(lambda (ind)
       (new-individual (funcall chr-mutator (individual-chromosome ind))
		       chr-evaluator))))
;;
;; This can be a simple template for further customization.
(defparameter *EC-output* *standard-output*)

(defun generational-elitism-EC (&key (population-size 500)
				chr-init chr-evaluator
				(chr-crossoveror nil) (chr-mutator nil)
				(ind-crossoveror nil) (ind-mutator nil)
				chr-printer
				preparer selector (better #'higher-fitness-than)
				(p-crossover 0.9) (p-mutation 0.05) (generations 50) (report-fitness-record nil))
  ;; if provided, chr-crossoveror is chromosome crossover binary function. Must be provided if ind-crossoveror is not provided.
  ;; if provided, chr-mutator is chromosome mutation unary function. Must be provided if ind-mutator is not provided.
  ;; if provided, ind-crossoveror is individual crossover binary function. Must be provided if chr-crossoveror is not provided.
  ;; if provided, ind-mutator is individual mutation unary function. Must be provided if chr-mutator is not provided.
  ;; returns the best found, and the last generation
  (let* ((n-crossover (ceiling (* population-size p-crossover)))
	 (n-mutation (ceiling (* population-size p-mutation)))
	 (n-reproduction (max 0 (- population-size n-crossover n-mutation 1)))
	 (cur-pop (fill-all-with (make-array population-size)
		    (new-individual (funcall chr-init) chr-evaluator)))
	 (best-so-far (best-of cur-pop better))
	 (next-pop (make-array population-size :initial-element nil)))
    ;;
    (or ind-crossoveror (setf ind-crossoveror #'(lambda (ind1 ind2)
						  (new-individual (funcall chr-crossoveror
									   (individual-chromosome ind1)
									   (individual-chromosome ind2))
								  chr-evaluator))))
    (or ind-mutator (setf ind-mutator #'(lambda (ind)
					  (new-individual (funcall chr-mutator
								   (individual-chromosome ind))
							  chr-evaluator))))
    ;;
    (dotimes (g generations)
      ;; next generation
      (let ((sel (funcall preparer cur-pop))
	    (cur-best (best-of cur-pop better)))
	(fill-with next-pop
	  (1 cur-best)
	  (n-crossover (funcall ind-crossoveror (funcall selector sel) (funcall selector sel)))
	  (n-mutation (funcall ind-mutator (funcall selector sel)))
	  (n-reproduction (funcall selector sel)))
	;; sort by fitness
	(sort next-pop better)
	;; report
	(if (funcall better cur-best best-so-far) (setf best-so-far cur-best))
	(format *EC-output* "===Generation ~d ===~%" (1+ g))
	(format *EC-output* "~%best-so-far: fitness: ~a " (individual-fitness best-so-far))
	(funcall chr-printer (individual-chromosome best-so-far))
	(terpri *EC-output*)
	;; fitness record
	(when report-fitness-record
	  (report-fitness-changes (1+ g))))
      ;;
      (let ((tmp cur-pop))
	(setf cur-pop next-pop)
	(setf next-pop tmp)))
    ;;
    (values best-so-far cur-pop)))
;;;
;;
