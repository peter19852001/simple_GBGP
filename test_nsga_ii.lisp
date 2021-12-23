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
;; But follow and adapt the code from https://gist.github.com/Tiagoperes/1779d5f1c89bae0cfdb87b1960bba36d
;;
;; found from this post: https://stackoverflow.com/a/40087950

(defparameter *eta-c* 1.0)

(defun sbx-one-var (x1 x2 xL xU)
  "Simulated Crossover of two real values X1 and X2, and the lower bound is xL, upper bound is xU.
*ETA-C* is the distribution index for the crossover
Returns the two children values."
  (if (<= (abs (- x1 x2)) 1.0e-14)
      ;; if the two number are the same, no crossover needed
      (values x1 x2)
      (let* ((y1 (min x1 x2))
             (y2 (max x1 x2))
             (r (random 1.0)))
        (flet ((one-child (beta)
                 (let* ((alpha (- 2.0 (expt beta (- (+ 1.0 *eta-c*)))))
                        (betaq (expt (if (<= r (/ 1.0 alpha))
                                         (* r alpha)
                                         (/ 1.0 (- 2.0 (* r alpha))))
                                     (/ 1.0 (+ 1.0 *eta-c*)))))
                   (* 0.5 (- (+ y1 y2)
                             (* betaq (- y2 y1))))
                   )))
          (values
           ;; two different beta for the two children, otherwise they
           ;; are the same, and share the same r.
           (one-child (+ 1.0 (/ (* 2.0 (- y1 xL))
                                (- y2 y1))))
           (one-child (+ 1.0 (/ (* 2.0 (- xU y2))
                                (- y2 y1))))))
      )))

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
             (+ 1 (/ (* 9 (loop :for i :from 1 :below n
                               :sum (aref x i)))
                     (- n 1))))))
    (vector (aref chr 0)
            (let ((gx (g chr)))
              (* gx (- 1 (expt (/ (aref chr 0)
                                  gx)
                               2)))))
    ))

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
