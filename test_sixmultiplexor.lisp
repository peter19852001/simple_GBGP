;;; 
;;; To test the Lisp version of GBGP on SixMultiplexor, using the PAGE strategy:
;;;   do 50 runs, if can find optimum in at least 90% of the runs, give the
;;;   mean and sd of the number of evaluations to reach optimum
;;;

;;; should load simple_ec.lisp, GBGP.lisp and sixmultiplexor.lisp first
(in-package :test-6-multiplexor)

(defun one-run (&key pop-size depth)
  ;; return number of evaluations to optimum if found.
  ;; return nil if optimum not found
  (format t ";;; Test pop-size: ~a	depth: ~a~%" pop-size depth)
  (multiple-value-bind (best n-evals-to-optimum)
      (solve-6-multiplexor :pop-size pop-size :depth depth)
    n-evals-to-optimum))

(defun mean-sd (lst)
  (if (null lst)
      (values "nan" "nan")
      (let ((m 0)
	    (m2 0)
	    (n 0))
	(dolist (x lst)
	  (incf n)
	  (incf m x)
	  (incf m2 (* x x)))
	(let* ((mean (/ m n))
	       (sd (sqrt (- (/ m2 n) (* mean mean)))))
	  (values (coerce mean 'float) (coerce sd 'float))))))

(defmacro open-file-no-abort ((stream filespec &rest options) &body body)
  ;; to prevent the accumulated output being deleted by close if error occurs, as in with-open-file
  `(let ((,stream (open ,filespec ,@options)))
     (unwind-protect
	  (progn ,@body)
       (when ,stream (close ,stream)))))


;;; changed not to use tree size in individual selection
(defun test ()
  (open-file-no-abort (*standard-output* "test_sixmultiplexor_dump.txt" :direction :output :if-exists :supersede)
    (open-file-no-abort (all-res "test_sixmultiplexor_all_res.txt" :direction :output :if-exists :supersede)
      (open-file-no-abort (summary "test_sixmultiplexor_summary.txt" :direction :output :if-exists :supersede)
	(format all-res "depth	popsize	run	nevals~%")
	(format summary "depth	pop_size	nsuccess	avg_eval	sd_eval	avg_best90	sd_best90	avg_worst90	sd_worst90~%")
	(dolist (pop-size '(50 250 500 1000 2000 3000))
	  (dolist (depth '(3 4 5 6 7 8 9 10))
	    (let ((lst-n-evals-to-optimum nil))
	      (dotimes (i 50)
		(let ((res (one-run :pop-size pop-size :depth depth)))
		  (if res (push res lst-n-evals-to-optimum))
		  (format all-res "~a	~a	~a	~a~%"
			  depth pop-size (1+ i) (if res res "fail"))))
	      ;;
	      (let* ((sorted-lst (sort (copy-list lst-n-evals-to-optimum) #'<))
		     (n-success (length sorted-lst))
		     (best45 (if (> n-success 45)
				 (subseq sorted-lst 0 45)
				 sorted-lst))
		     (worst45 (if (> n-success 45)
				  (nthcdr (- n-success 45) sorted-lst)
				  sorted-lst)))
		(multiple-value-bind (avg sd) (mean-sd sorted-lst)
		  (multiple-value-bind (avg-best45 sd-best45) (mean-sd best45)
		    (multiple-value-bind (avg-worst45 sd-worst45) (mean-sd worst45)
		      (format summary "~a	~a	~a	~a	~a	~a	~a	~a	~a~%"
			      depth pop-size n-success avg sd
			      avg-best45 sd-best45 avg-worst45 sd-worst45))))))))))))
;;;
;; (test)
;;;;

