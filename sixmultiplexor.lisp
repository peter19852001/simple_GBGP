;;;
;;; to test the Grammar Based Genetic Programming framework
;;;
;;; test on the simple 6-Multiplexor problem, using grammar described in
;;;  P.A. Whigham, Grammatically-based genetic programming.
;;;    in "Proceedings of the Workshop on Genetic Programming: From Theory to Real-World Applications,
;;;    ed. by J.P. Rosca (Tahoe City, California, USA, 1995), pp. 33--41.
;;; Briefly, the 6-Multiplexor is:
;;;   Given boolean inputs a0, a1, d0, d1, d2, d3,
;;;   Use operators (IF X Y Z), (AND X Y), (OR X Y), (NOT X), to output
;;;     d0 if a0=0, a1=0;
;;;     d1 if a0=0, a1=1;
;;;     d2 if a0=1, a1=0;
;;;     d3 if a0=1, a1=1;

;; should load simple_ec.lisp, and GBGP.lisp first

(defpackage :test-6-multiplexor
  (:use :simple-grammar-based-GP :simple-evolutionary-computation :common-lisp))
(in-package :test-6-multiplexor)

(defparameter *6-multiplexor-grammar-external*
  ;; <S> is the starting symbol
  `((<S> :> <B>)
    (<B> -> and <B> <B>
	 :or or <B> <B>
	 :or not <B>
	 :or if <B> <B> <B>)
    (<B> :> <T>)
    (<T> :> a0 :or a1 :or d0 :or d1 :or d2 :or d3)))

(defparameter *6-multiplexor-grammar* (grammar-to-internal *6-multiplexor-grammar-external*))

(defparameter *init-max-tree-height* 6)
(defun random-6-multiplexor () (limit-expand-non-terminal-from-grammar *6-multiplexor-grammar* '<S> *init-max-tree-height*))

(defun correct-6-multiplexor (a0 a1 d0 d1 d2 d3)
  (if a0
      (if a1 d3 d2)
      (if a1 d1 d0)))

(defun compile-6-multiplexor-chr (chr)
  "Turn the parse-tree into Lisp S-expression and eval it to give a function."
  (eval `#'(lambda (a0 a1 d0 d1 d2 d3)
	     (declare (ignorable a0 a1 d0 d1 d2 d3))
	     ,(genotype-to-phenotype chr))))

(defparameter *n-evals-to-optimum* nil)
(defparameter *n-evals* 0)
(defun count-correct-6-multiplexor (chr)
  (incf *n-evals*)
  (let ((chr-f (compile-6-multiplexor-chr chr))
	(correct-count 0))
    (dolist (a0 '(t nil))
      (dolist (a1 '(t nil))
	(dolist (d0 '(t nil))
	  (dolist (d1 '(t nil))
	    (dolist (d2 '(t nil))
	      (dolist (d3 '(t nil))
		(if (eq (funcall chr-f a0 a1 d0 d1 d2 d3) (correct-6-multiplexor a0 a1 d0 d1 d2 d3))
		    (incf correct-count))))))))
    (when (= correct-count 64) (setf *n-evals-to-optimum* *n-evals*))
    correct-count))

(defun solve-6-multiplexor (&key (pop-size 500) (depth 6))
  ;; return the best phenotype, number of evaluations to reach optimum (if reached, nil otherwise)
  (let ((*n-evals* 0)
	(*n-evals-to-optimum* nil))
    (multiple-value-bind (best best-ind n-evals last-pop)
	(steady-state-elitism-plain-GBGP *6-multiplexor-grammar* '<S>  #'count-correct-6-multiplexor
					 :population-size pop-size :p-crossover 0.7 :p-mutation 0.2 :generations 20000 :init-max-height depth :max-tree-height depth
					 :max-evals 200000 :optimum 64 :stop-at-optimum t)
      (values best *n-evals-to-optimum*))))

;(solve-6-multiplexor)
;; at least seems to work, though the obtained function has too many nodes. Probably the tree size need to be controlled better.
