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
        (cond ((eql v1 v2) nil)
              ((funcall (aref better i) v1 v2)
               (setf has-strictly-better))
              (t (return nil)))))))

(defun individual-dominator (better)
  #'(lambda (ind-dominating ind-dominated)
      (fitness-dominate-p
       (individual-fitness ind-dominating)
       (individual-fitness ind-dominated)
       better)))

;;; 
(defun fast-non-dominated-sort (pop ind-dominate-p)
  "Calculate the list of list of non-dominated MO-INDIVIDUAL with increasing ranks.
POP is a vector of INDIVIDUAL.
IND-DOMINATE-P is (lambda (ind-dominating ind-dominated) ...)."
  (let ((mo-pop (map 'list
                     #'(lambda (x) (make-mo-ind :ind x))
                     pop))
        (n-inds (length pop))
        (first-front nil)
        (all-fronts nil))
    (dolist (p mo-pop)      
      (dolist (q mo-pop)
        (when (funcall ind-dominate-p (mod-ind-ind p) (mo-ind-ind q))
          ;; p dominates q
          (push q (mo-ind-dominated-by p))
          (incf (mo-ind-n-dominating q)))))
    (dolist (p mo-pop)
      (when (= 0 (mo-ind-n-dominating p))
        ;; p belongs to the first front
        (setf (mo-ind-rank p) 1)
        (push p first-front)))
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
      (push cur-front all-fronts))
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
        (loop :for i :from 2 :below n
             :do (incf (mo-ind-distance (aref sort-by-obj i))
                       (/ (- (mo-fitness-obj-i (aref sort-by-obj (1+ i)) obj)
                             (mo-fitness-obj-i (aref sort-by-obj (- i 1)) obj))
                          obj-value-range)))))))

(defun NSGA-II-next-gen-parents (pop ind-dominate-p n)
  "POP is a vector of (more than N) INDIVIDUAL, select N INDIVIDUAL
according to nondomination and crowding distance."
  (let ((new-parents (make-array n))
        (n-inds 0))
    (do ((non-dom-inds (fast-non-dominated-sort pop ind-dominate-p) (cdr non-dom-inds)))
        ((> (+ n-inds (length (car non-dom-inds))) n)
         ;; the current front will fill-up, so take only as many as needed
         ;; all the mo-inds in the front have the same rank, so sort by decreasing crowding distance
         (do ((last-front (sort (crowding-distance-assignment (car non-dom-inds))
                                 #'> :key #'mo-ind-distance)
                          (cdr last-front)))
             ((>= n-inds n))
           (setf (aref new-parents n-inds) (car last-front))
           (incf n-inds)))
      ;; still not filled up yet, calculate crowding distance for later use
      (crowding-distance-assignment (car non-dom-inds))
      ;; put the mo-inds as new parents
      (dolist (x (car non-dom-inds))
        (setf (aref new-parents n-inds) x)
        (incf n-inds)))
    ;;
    new-parents))

(defun NSGA-II ()
  ;; TODO
  )
