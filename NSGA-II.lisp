;;; 
;;; A simple implementation of NSGA-II without constraints, following
;;; 
;;; Deb, K., Pratap, A., Agarwal, S., & Meyarivan,
;;; T. A. M. T. (2002). A fast and elitist multiobjective genetic
;;; algorithm: NSGA-II. IEEE transactions on evolutionary computation,
;;; 6(2), 182-197.
;;;
;;; https://www.cse.unr.edu/~sushil/class/gas/papers/nsga2.pdf

(defpackage :simple-nsga-II
  (:use :common-lisp :simple-evolutionary-computation)
  (:nicknames :nsga2)
  (:export
   :NSGA-II
   :MO-IND-IND
   :output-pareto-front-fitness-as-csv
   :output-pareto-front-fitness-to-csv-file)
  )

(in-package :simple-nsga-II)

;;; use the same individual as :EC, but add a wrapper for extra tempoary attribute
(defstruct mo-ind
  "Wraps an INDIVIDUAL, and add extra attributes for some calculations"
  ind ;; the wrapped INDIVIDUAL
  (n-dominating 0) ;; number of individuals dominating this one
  (dominated-by nil) ;; list of MO-IND dominated by this one
  (rank nil)
  ;; for crowding distance
  distance
  )

;; convenience functions
(defun mo-fitness (mo)
  (individual-fitness (mo-ind-ind mo)))

(defun mo-fitness-obj-i (mo i)
  "i-th objective of MO-IND mo"
  (aref (mo-fitness mo) i))

;;; fitness domination
(defun fitness-dominate-p (fitness-dominating fitness-dominated better)
  "Whether FITNESS-DOMINATING dominates FITNESS-DOMINATED, where the values are compared by the operators in BETTER.
All FITNESS-DOMINATING, FITNESS-DOMINATED and BETTER are vectors."
  (let ((has-strictly-better nil))
    (dotimes (i (length better) has-strictly-better)
      (let ((v1 (aref fitness-dominating i))
            (v2 (aref fitness-dominated i)))
        (cond ((or (not (numberp v1))
                   (not (numberp v2)))
               (return nil))
              ((eql v1 v2) nil)
              ((funcall (aref better i) v1 v2)
               (setf has-strictly-better t))
              (t (return nil)))))))

(defun individual-dominator (better)
  #'(lambda (ind-dominating ind-dominated)
      (fitness-dominate-p
       (individual-fitness ind-dominating)
       (individual-fitness ind-dominated)
       better)))

(defun mo-ind-better-p (better-ind worse-ind)
  ;; prefer either lower rank, or if same rank, larger crowding distance
  (or (< (mo-ind-rank better-ind)
         (mo-ind-rank worse-ind))
      (and (= (mo-ind-rank better-ind)
              (mo-ind-rank worse-ind))
           (> (mo-ind-distance better-ind)
              (mo-ind-distance worse-ind)))))

(defun population-summary (pop)
  "POP is a vector of MO-INDIVIDUALs, print their fitness and chrosome."
  (dotimes (i (length pop))
    (format t "~A: ~A~%"
            i (mo-ind-ind (aref pop i)))))
;;; 
(defun fast-non-dominated-sort (pop ind-dominate-p)
  "Calculate the list of list of non-dominated MO-INDIVIDUAL with increasing ranks.
POP is a vector of MO-INDIVIDUAL and their ranks will be updated.
The rank of the first front is 1.
IND-DOMINATE-P is (lambda (ind-dominating ind-dominated) ...)."
  ;;(format t "At fast-non-dominated-sort~%")
  ;;(population-summary pop)
  (let ((n-inds (length pop))
        (first-front nil)
        (all-fronts nil))
    ;; first clear the auxiliary info
    (dotimes (i n-inds)
      (let ((p (aref pop i)))
        (setf (mo-ind-n-dominating p) 0
              (mo-ind-dominated-by p) nil
              (mo-ind-rank p) nil)))
    ;;
    (dotimes (i n-inds)
      (dotimes (j n-inds)
        (let ((p (aref pop i))
              (q (aref pop j)))
          (when (funcall ind-dominate-p (mo-ind-ind p) (mo-ind-ind q))
            ;; p dominates q
            (push q (mo-ind-dominated-by p))
            (incf (mo-ind-n-dominating q))))))
    ;;
    (dotimes (i n-inds)
      (let ((p (aref pop i)))
        (when (= 0 (mo-ind-n-dominating p))
          ;; p belongs to the first front
          (setf (mo-ind-rank p) 1)
          (push p first-front))))
    (push first-front all-fronts)
    ;; other fronts
    (do ((cur-front nil nil)
         (prev-front first-front cur-front)
         (cur-rank 2 (1+ cur-rank)))
        ((null prev-front))
      (dolist (p prev-front)
        (dolist (q (mo-ind-dominated-by p))
          (when (= 0 (decf (mo-ind-n-dominating q)))
            (setf (mo-ind-rank q) cur-rank)
            (push q cur-front))))
      ;;
      (when cur-front
        (push cur-front all-fronts)))
    ;;
    (nreverse all-fronts)))

(defun crowding-distance-assignment (inds)
  "Calculate crowding distance assignment and will modify the
MO-INDIVIDUALs in INDS which is a list of MO-INDIVIDUAL"
  (let ((v-inds (coerce inds 'vector))
        (n (length inds))
        ;; number of objectives, use the fitness of first individual to determine
        (n-objs (length (mo-fitness (first inds)))))
    (dolist (i inds) (setf (mo-ind-distance i) 0))
    (dotimes (obj n-objs)
      (let* ((sort-by-obj
              (sort v-inds #'< :key #'(lambda (x) (mo-fitness-obj-i x obj))))
             (smallest-ind (aref sort-by-obj 0))
             (largest-ind (aref sort-by-obj (- n 1)))
             (obj-value-range (+ 1e-15
                                 (- (mo-fitness-obj-i largest-ind obj)
                                    (mo-fitness-obj-i smallest-ind obj)))))
        ;; sbcl specific, for positive infinity
        (setf (mo-ind-distance smallest-ind) SB-EXT:DOUBLE-FLOAT-POSITIVE-INFINITY
              (mo-ind-distance largest-ind) SB-EXT:DOUBLE-FLOAT-POSITIVE-INFINITY)
        ;;
        (loop :for i :from 1 :below (- n 1)
             :do (incf (mo-ind-distance (aref sort-by-obj i))
                       (/ (- (mo-fitness-obj-i (aref sort-by-obj (1+ i)) obj)
                             (mo-fitness-obj-i (aref sort-by-obj (- i 1)) obj))
                          obj-value-range)))))))

(defun NSGA-II-next-gen-parents (pop ind-dominate-p n &optional cur-gen report-func)
  "POP is a vector of (more than N) INDIVIDUAL, select N INDIVIDUAL
according to nondomination and crowding distance.
REPORT-FUNC is the same as in NSGA-II to report the current non-dominated front for generation CUR-GEN."
  ;;(format t "At NSGA-II-next-gen-parents~%")
  ;;(population-summary pop)
  (let ((new-parents (make-array n))
        (n-inds 0)
        (all-fronts (fast-non-dominated-sort pop ind-dominate-p)))
    (when report-func
      (funcall report-func cur-gen (car all-fronts)))
    ;;
    (do ((non-dom-inds all-fronts (cdr non-dom-inds)))
        ((> (+ n-inds (length (car non-dom-inds))) n)
         ;; the current front will fill-up, so take only as many as needed
         ;; all the mo-inds in the front have the same rank, so sort by decreasing crowding distance
         (format t "Caclulate crowding distances of last front of ~A inds.~%" (length (car non-dom-inds)))
         (finish-output)
         (crowding-distance-assignment (car non-dom-inds))
         (format t "~A inds already filled, get the last front.~%" n-inds)
         (do ((last-front (sort (car non-dom-inds)
                                 #'> :key #'mo-ind-distance)
                          (cdr last-front)))
             ((>= n-inds n))
           (setf (aref new-parents n-inds) (car last-front))
           (incf n-inds)))
      ;; still not filled up yet, calculate crowding distance for later use
      (format t "Caclulate crowding distances of ~A inds.~%" (length (car non-dom-inds)))
      (finish-output)
      (crowding-distance-assignment (car non-dom-inds))
      ;; put the mo-inds as new parents
      (dolist (x (car non-dom-inds))
        (setf (aref new-parents n-inds) x)
        (incf n-inds)))
    ;;
    new-parents))

(defun tournament-select-one (pop)
  "POP is a vector of MO-INDIVIDUALs already with rank and crowding
  distance calculated, select one as potential parent through
  tournament."
  (let* ((n (length pop))
         (ind-i (aref pop (random n)))
         (ind-j (aref pop (random n))))
    (if (mo-ind-better-p ind-i ind-j)
        ind-i
        ind-j)))

(defun NSGA-II (&key (population-size 500)
				  chr-init chr-evaluator
				  (chr-crossoveror nil) (chr-mutator nil)
				  ;; chr-printer
                  better 
				  (p-mutation 0.05)
                  (generations 50)
                  (report-func nil))
  "The simple version of NSGA-II without constraint.

POPULATION-SIZE: population size.
CHR-INIT: nullary function to generate a chromosome.
CHR-EVALUATOR: unary function to evaluate a chromosome to get the fitness which is a vector of values for the objectives.
CHR-CROSSOVEROR: binary function to crossover two chromosomes and return a pair of chromosomes.
CHR-MUTATOR: unary function to mutate a chromosome and return a mutated chromosome.
CHR-PRINTER: unary function to print a chromosome in a readable way.
BETTER: a vector of binary predicate to indicate what is better for each of the multiple objectives; e.g. #'< indicates an objective should be minimized.
P-MUTATION: probability of mutation.
GENERATIONS: number of generations.
REPORT-FUNC: if non-nil, should be a function of (lambda (generation mo-inds) ...) that reports the current front as appropriate.

Refer to 
Deb, K., Pratap, A., Agarwal, S., & Meyarivan, T. A. M. T. (2002). A fast and elitist multiobjective genetic algorithm: NSGA-II. IEEE transactions on evolutionary computation, 6(2), 182-197.
"
  (labels ((gen-offsprings (pop)
             ;; pop is a vector of MO-INDIVIDUALs with rank and crowding distance
             ;; to produce a vector of the same length of MO-INDIVIDUALs as offsprings
             (let* ((n (length pop))
                    (out (make-array n)))
               (do ((i 0 (+ i 2)))
                   ((>= i n) out)
                 (let* ((p1 (tournament-select-one pop))
                        (p2 (tournament-select-one pop))
                        ;; chr-crossoveror should produce a pair of two chromosomes
                        (children-chrs
                         (funcall chr-crossoveror
                                  (individual-chromosome (mo-ind-ind p1))
                                  (individual-chromosome (mo-ind-ind p2))))
                        (chr1 (if (bernoulli p-mutation)
                                  (funcall chr-mutator (car children-chrs))
                                  (car children-chrs)))
                        (chr2 (if (bernoulli p-mutation)
                                  (funcall chr-mutator (cdr children-chrs))
                                  (cdr children-chrs)))
                        (ind1 (new-individual chr1 chr-evaluator))
                        (ind2 (new-individual chr2 chr-evaluator)))
                   (setf (aref out i)
                         (make-mo-ind :ind ind1))
                   (when (< (1+ i) n)
                     (setf (aref out (1+ i))
                           (make-mo-ind :ind ind2)))
                   )))))
    ;;
    (let* ((ind-dominate-p (individual-dominator better))
           (cur-parents
            ;; a vector of MO-INDIVIDUALs
            (progn
              (format t "Generate initial population.~%")
              (finish-output)
              (fill-all-with (make-array population-size)
                             (make-mo-ind :ind (new-individual
                                                (funcall chr-init)
                                                chr-evaluator)))))
           (non-dom-fronts
            (progn
              (format t "Calculate initial non-dominated pareto fronts.~%")
              (finish-output)
              (fast-non-dominated-sort cur-parents ind-dominate-p))))
      ;; initial population's handling is a little different
      ;; first get ranks and crowding distance for initial parents
      (format t "Calculate initial crowding distances.~%")
      (finish-output)
      (dolist (inds non-dom-fronts)
        (crowding-distance-assignment inds))
      ;;
      (dotimes (g generations)
        ;; next generation
        (format t "Generation ~A~%Generate offsprings~%" g)
        (finish-output)
        ;; first generate offsprings
        (let* ((offsprings (gen-offsprings cur-parents))
               (new-pop (concatenate 'vector
                                     cur-parents offsprings)))
          (format t "Determine parents of next generation.~%")
          (finish-output)
          (setf cur-parents (NSGA-II-next-gen-parents
                             new-pop
                             ind-dominate-p
                             population-size
                             g report-func))))
      ;; return the non-dominated solutions.
      (format t "Return the final front.~%")
      ;; now cur-parents has ranks but we have lost the fronts because
      ;; we did not keep track of them, so simply loop through it to
      ;; get the non-dominated set.
      (let ((final-non-dom nil))
        (dotimes (i (length cur-parents))
          (let ((p (aref cur-parents i)))
            (when (= 1 (mo-ind-rank p))
              (push p final-non-dom))))
        ;;
        (when report-func
          (funcall report-func generations final-non-dom))
        ;;
        final-non-dom)
      )))

(defun output-pareto-front-fitness-as-csv (inds)
  (format t "individual")
  (let ((n-objs (length (mo-fitness (first inds)))))
    (dotimes (j n-objs)
      (format t ",objective~A" j))
    (terpri)
    ;;
    (flet ((print-fitnesses (fitness)
             (dotimes (j (length fitness))
               (format t ",~A" (aref fitness j)))))
      (do ((i 0 (1+ i))
           (ds inds (cdr ds)))
          ((null ds))
        (format t "~A" i)
        (print-fitnesses (mo-fitness (car ds)))
        (terpri)))))

(defun output-pareto-front-fitness-to-csv-file (inds out-file)
  (with-open-file (*standard-output* out-file
                                     :direction :output
                                     :if-exists :supersede)
    (output-pareto-front-fitness-as-csv inds)))
