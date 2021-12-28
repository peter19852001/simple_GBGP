;;;
;;; To test NSGA-II on selected problems mentioned in the paper.
;;;

;;; should load simple_ec.lisp, NSGA-II.lisp.
(in-package :simple-nsga-ii)

;; we use a real number or a vector of real numbers for the
;; chromosome. Following the NSGA-II paper, we use simulated binary
;; crossover and polynomial mutation as described in
;;
;; Deb, Kalyanmoy, and Ram Bhushan Agrawal. "Simulated binary
;; crossover for continuous search space." Complex systems 9.2 (1995):
;; 115-148.
;;
;; But follow and adapt the code from https://gist.github.com/Tiagoperes/1779d5f1c89bae0cfdb87b1960bba36d for the bound version.
;;
;; found from this post: https://stackoverflow.com/a/40087950
;;
;; But we do not use the bounded version here.
;;
;;
;; For unbound version, follows the formulae in the paper:
;; 
;; Deb, Kalyanmoy, and Hans-Georg Beyer. "Self-adaptive genetic
;; algorithms with simulated binary crossover." Evolutionary
;; computation 9.2 (2001): 197-221.

;; For mutation, NSGA-II paper mentions polynomial operator, but seems to give the wrong citation, so search and found mentioning of polynomial mutation operator for real-coded GA in:
;;
;; Deb, Kalyanmoy, and Debayan Deb. "Analysing mutation schemes for
;; real-parameter genetic algorithms." International Journal of
;; Artificial Intelligence and Soft Computing 4.1 (2014): 1-28.
;;
;; which listed the formula proposed in an earlier paper.

(defparameter *eta-c* 20.0
  "Distribution index of Simulated Binary Crossover operator used in NSGA-II paper.")
(defparameter *eta-m* 20.0
  "Distribution index of Simulated Binary Mutation operator used in NSGA-II paper.")

(defun sbx-one-var-b (x1 x2 xL xU eta-c)
  "Simulated Crossover of two real values X1 and X2, and the lower bound is XL, upper bound is XU.
ETA-C is the distribution index for the crossover.
Returns the two children values."
  (if (<= (abs (- x1 x2)) 1.0e-14)
      ;; if the two number are the same, no crossover needed
      (values x1 x2)
      (let* ((y1 (min x1 x2))
             (y2 (max x1 x2))
             (r (random 1.0)))
        (flet ((one-child (beta)
                 (let* ((alpha (- 2.0 (expt beta (- (+ 1.0 eta-c)))))
                        (betaq (expt (if (<= r (/ 1.0 alpha))
                                         (* r alpha)
                                         (/ 1.0 (- 2.0 (* r alpha))))
                                     (/ 1.0 (+ 1.0 eta-c)))))
                   (* 0.5 (- (+ y1 y2)
                             (* betaq (- y2 y1))))
                   ))
               (clip-to-range (z)
                 (cond ((< z xL) xL)
                       ((> z xU) xU)
                       (t z))))
          (values
           ;; two different beta for the two children, otherwise they
           ;; are the same, and share the same r.
           (clip-to-range
            (one-child (+ 1.0 (/ (* 2.0 (- y1 xL))
                                 (- y2 y1)))))
           (clip-to-range
            (one-child (+ 1.0 (/ (* 2.0 (- xU y2))
                                 (- y2 y1)))))))
        )))

(defun sbx-one-var-ub (x1 x2 eta)
  "Simulated Binary Crossover for real values X1 and X2 without bounds.
ETA is the distribution index.
Returns two values as children."
  (let* ((u (random 1.0))
         (betaq (expt
                 (if (<= u 0.5)
                     (* 2.0 u)
                     (/ 1.0 (* 2.0 (- 1.0 u))))
                 (/ 1.0 (+ 1.0 eta))))
         (a (+ 1.0 betaq))
         (b (- 1.0 betaq)))
    (values
     ;; child 1
     (* 0.5 (+ (* a x1) (* b x2)))
     ;; child 2
     (* 0.5 (+ (* b x1) (* a x2))))))

(defun sbx-vars-b (x1 x2 xL xU eta-c)
  "Simulated Crossover of two real vectors X1 and X2 of the same length, and the lower bound vector is XL, upper bound vector is XU.
ETA-C is the distribution index for the crossover.
Returns the two children values as a pair of vectors."
  (let* ((n (length x1))
         (c1 (make-array n))
         (c2 (make-array n)))
    (dotimes (i n)
      (let ((p1 (aref x1 i))
            (p2 (aref x2 i)))
        (if (<= (random 1.0) 0.5)
            ;; crossover
            (multiple-value-bind (y1 y2)
                (sbx-one-var-b p1 p2 (aref xL i) (aref xU i) eta-c)
              (setf (aref c1 i) y1
                    (aref c2 i) y2))
            ;; just copy from parent
            (setf (aref c1 i) p1
                  (aref c2 i) p2))))
    (cons c1 c2)))

(defun sbx-vars-ub (x1 x2 eta-c)
  "Simulated Crossover of two real vectors X1 and X2 of the same length, the unbounded version.
ETA-C is the distribution index for the crossover.
Returns the two children values as a pair of vectors."
  (let* ((n (length x1))
         (c1 (make-array n))
         (c2 (make-array n)))
    (dotimes (i n)
      (let ((p1 (aref x1 i))
            (p2 (aref x2 i)))
        (if (<= (random 1.0) 0.5)
            ;; crossover
            (multiple-value-bind (y1 y2)
                (sbx-one-var-ub p1 p2 eta-c)
              (setf (aref c1 i) y1
                    (aref c2 i) y2))
            ;; just copy from parent
            (setf (aref c1 i) p1
                  (aref c2 i) p2))))
    (cons c1 c2)))

(defun poly-mutate-one-var (p xL xU eta)
  "Polynomial Mutation Operator for Real-Coded GA, for real value P, where XL and XU are the lower and upper bounds respectively. ETA is the distribution index."
  (let* ((u (random 1.0))
         (child (if (<= u 0.5)
                    ;;
                    (+ p (* (- (expt (* 2.0 u)
                                     (/ 1.0 (+ 1.0 eta)))
                               1.0)
                            (- p xL)))
                    ;;
                    (+ p (* (- 1.0
                               (expt (* 2.0 (- 1.0 u))
                                     (/ 1.0 (+ 1.0 eta))))
                            (- xU p))))))
    ;; clip to range
    (cond ((< child xL) xL)
          ((> child xU) xU)
          (t child))))

(defun poly-mutate-vars (x xL xU eta &optional pm)
  "Polynomial Mutation Operator for Real-Coded GA applied to real vector X, where XL and XI are the lower and upper bounds respectively. ETA is the distribution index.
Use a simple schemme of mutating each entry of X independently with a probability of PM."
  ;; use explicit loop for simplicity
  (let* ((n (length x))
         (v (make-array n)))
    (if (null pm) (setf pm (/ 1.0 n)))
    (dotimes (i n v)
      (setf (aref v i)
            (if (<= (random 1.0) pm)
                ;; mutate
                (poly-mutate-one-var
                 (aref x i)
                 (aref xL i)
                 (aref xU i)
                 eta)
                ;; just copy
                (aref x i))))))

;;;;

;;;;;;;;
;; SCH
;; chromosome is one real value in [-10^3, 10^3]
(defun sch-chr-init ()
  (- (random 2000.0) 1000.0))

(defun sch-fitness (chr)
  "chr is just one real value x.
Two objective functions:
f1(x) = x^2
f2(x) = (x-2)^2"
  (vector (expt chr 2)
          (expt (- chr 2) 2)))

(defun sch-chr-crossover (chr1 chr2)
  "CHR1 and CHR2 are real values."
  (multiple-value-bind (c1 c2)
      (sbx-one-var-ub chr1 chr2
                      *eta-c*)
    (cons c1 c2)))

(defun sch-chr-mutate (chr)
  (poly-mutate-one-var chr
                       -1000.0 1000.0
                       *eta-m*))

(defparameter *sch-better* (vector #'< #'<))

(defparameter *sch-fronts*
  (NSGA-II :population-size 100
           :chr-init #'sch-chr-init
           :chr-evaluator #'sch-fitness
           :chr-crossoveror #'sch-chr-crossover
           :chr-mutator #'sch-chr-mutate
           :better *sch-better*
           :p-mutation 1.0
           :generations 250))

(output-pareto-front-fitness-to-csv-file *sch-fronts* "sch-test-fronts.csv")
;; the plot of fitness (using R) looks similar to that shown in the paper.

;;;;;;;;
;; KUR
;; chromosome is a real vector of length 3 in [-5, 5]
(defun kur-chr-init ()
  (vector (- (random 10.0) 5.0)
          (- (random 10.0) 5.0)
          (- (random 10.0) 5.0)))

(defun kur-fitness (chr)
  "chr is a vector of length n.
f1(x) = \sum_{i=1}^{n-1}(-10 exp(-0.2 \sqrt{x_i^2 + x_{i+1}^2}))
f2(x) = \sum_{i=1}^{n-1}(|x_i|^0.8 + 5sin{x_i^3}"
  (vector
   ;; f1
   (loop :for i :below (- (length chr) 1)
      :sum (* -10 (exp (* -0.2 (sqrt (+ (expt (aref chr i) 2)
                                        (expt (aref chr (1+ i)) 2)))))))
   ;; f2
   (loop :for i :below (length chr)
        :sum (+ (expt (abs (aref chr i)) 0.8)
                (* 5 (sin (expt (aref chr i) 3)))))))

(defparameter *kur-xL* #(-5.0 -5.0 -5.0))
(defparameter *kur-xU* #(5.0 5.0 5.0))

(defun kur-chr-crossover (chr1 chr2)
  "CHR1 and CHR2 are real vectors of length 3."
  (sbx-vars-ub chr1 chr2
               *eta-c*))

(defun kur-chr-mutate (chr)
  (poly-mutate-vars chr
                    *kur-xL* *kur-xU*
                    *eta-m*))

(defparameter *kur-better* (vector #'< #'<))

(defparameter *kur-fronts*
  (NSGA-II :population-size 100
           :chr-init #'kur-chr-init
           :chr-evaluator #'kur-fitness
           :chr-crossoveror #'kur-chr-crossover
           :chr-mutator #'kur-chr-mutate
           :better *kur-better*
           :p-mutation (/ 1.0 3.0)
           :generations 250))

(output-pareto-front-fitness-to-csv-file *kur-fronts* "kur-test-fronts.csv")
;; the plot of fitness (using R) looks similar to that shown in the paper.

;;;;;;;;
;; ZDT2
;; chromosome is a real vector of length 30 in [0,1]
(defun zdt2-chr-init ()
  (fill-all-with (make-array 30)
    (random 1.0)))

(defun zdt2-fitness (chr)
  "chr is a vector of length n.
f1(x) = x_1
f2(x) = g(x)[1 - (x_1 / g(x))^2]
g(x) = 1 + 9(\sum_{i=2}^n x_i)/(n-1)"
  (flet ((g (x)
           (let ((n (length x)))
             (+ 1.0 (/ (* 9.0 (loop :for i :from 1 :below n
                                 :sum (aref x i)))
                       (- n 1))))))
    (vector (aref chr 0)
            (let ((gx (g chr)))
              (* gx (- 1.0 (expt (/ (aref chr 0)
                                    gx)
                                 2)))))
    ))

(defparameter *zdt2-xL* (make-array 30 :initial-element 0.0))
(defparameter *zdt2-xU* (make-array 30 :initial-element 1.0))

(defun zdt2-chr-crossover (chr1 chr2)
  "CHR1 and CHR2 are real vectors of length 30."
  (sbx-vars-b chr1 chr2
              *zdt2-xL* *zdt2-xU*
              *eta-c*))

(defun zdt2-chr-mutate (chr)
  (poly-mutate-vars chr
                    *zdt2-xL* *zdt2-xU*
                    *eta-m*))

(defparameter *zdt2-better* (vector #'< #'<))

(defparameter *zdt2-fronts*
  (NSGA-II :population-size 100
           :chr-init #'zdt2-chr-init
           :chr-evaluator #'zdt2-fitness
           :chr-crossoveror #'zdt2-chr-crossover
           :chr-mutator #'zdt2-chr-mutate
           :better *zdt2-better*
           :p-mutation (/ 1.0 30.0)
           :generations 2500))

(output-pareto-front-fitness-to-csv-file *zdt2-fronts* "zdt2-test-fronts.csv")
;; the plot of fitness (using R) is not quite the same as the plot in
;; the paper, maybe the crossover and mutation operators are not
;; exactly the same.

;; seems zdt2 is binary-coded in the paper, that may explain the discrepancies of results.

;;;;;;;;
;; ZDT4
;; chromosome is a real vector of length 10, where x_1 is in [0,1], others are in [-5,5]
(defun zdt4-chr-init ()
  (fill-with (make-array 10)
    (1 (random 1.0))
    (9 (- (random 10.0) 5.0))))

(defun zdt4-fitness (chr)
  "chr is a vector of length n.
f1(x) = x_1
f2(x) = g(x)[1 - \sqrt{x_1 / g(x)}]
g(x) = 1 + 10(n-1) + \sum_{i=2}^n [x_i^2 - 10cos(4\pi x_i)]"
  (flet ((g (x)
           (let ((n (length x)))
             (+ 1
                (* 10 (- n 1))
                (loop :for i :from 1 :below n
                     :sum (- (expt (aref x i) 2)
                             (* 10 (cos (* 4 pi (aref x i))))))))
           ))
    (vector (aref chr 0)
            (let ((gx (g chr)))
              (* gx (- 1 (sqrt (/ (aref chr 0)
                                  gx))))))
    ))

(defparameter *zdt4-xL*
  (let ((v (make-array 10 :initial-element -5.0)))
    (setf (aref v 0) 0.0)
    v))
(defparameter *zdt4-xU*
  (let ((v (make-array 10 :initial-element 5.0)))
    (setf (aref v 0) 1.0)
    v))

(defun zdt4-chr-crossover (chr1 chr2)
  "CHR1 and CHR2 are real vectors of length 30."
  (sbx-vars-b chr1 chr2
              *zdt4-xL* *zdt4-xU*
              *eta-c*))

(defun zdt4-chr-mutate (chr)
  (poly-mutate-vars chr
                    *zdt4-xL* *zdt4-xU*
                    *eta-m*))

(defparameter *zdt4-better* (vector #'< #'<))

(defparameter *zdt4-fronts*
  (NSGA-II :population-size 100
           :chr-init #'zdt4-chr-init
           :chr-evaluator #'zdt4-fitness
           :chr-crossoveror #'zdt4-chr-crossover
           :chr-mutator #'zdt4-chr-mutate
           :better *zdt4-better*
           :p-mutation 0.1
           :generations 2500))

(output-pareto-front-fitness-to-csv-file *zdt4-fronts* "zdt4-test-fronts.csv")
;; the plot of fitness (using R) is not quite the same as that in the
;; paper, maybe because the crossover and mutation operators are
;; exactly the same.
