;;;
;;; A very simplistic application of Grammar Based Genetic Programming
;;; Useful reference:
;;; Robert I. Mckay, Nguyen Xuan Hoai, Peter Alexander Whigham, Yin Shan, and Michael O'Neill. 2010.
;;;  Grammar-based Genetic Programming: a survey. Genetic Programming and Evolvable Machines 11,
;;;  3-4 (September 2010), 365-396. DOI=10.1007/s10710-010-9109-y http://dx.doi.org/10.1007/s10710-010-9109-y 
;;;

;;; The external representation of a grammar is, simply represented as a list of lists,
;;; where each list represents a grammar rule.
;;; The basic form of a grammar rule is (N -> ...) where N (a symbol) represent non-terminals, and the second
;;; element should be the symbol ->, which is more for readability. ... are the terminals and non-terminals.
;;; The set of non-terminals is determined by the set of head symbol of the rules.
;;; At the end of a grammar rule, if the keyword :action appears, it should be followed by
;;;   a function of two variables: the rule (object) and the list of constructed objects for its non-terminal children.
;;;   This function should produce a value for the parse tree node for the current rule, when converting from genotype
;;;   to phenotype.
;;;   If action is not provided, the default is to #'collect-list, which collects the list of terminals and converted non-terminals of the rule.
;;; E.g. a rule (Sent -> Sub Pred)
;;; Short-form of a grammar rule:
;;;   (N :> ...) will make the defalt action #'collect-leaf, for which the rule body should contain only one term, and it will
;;;     be returned in rule action.
;;;   (N -> ... :or ... :or ...) or (N :> ... :or ... :or ...), :or can be used to separate alternatives to the non-terminal N,
;;;     and each alternative can have its own rule action. 
;;;   

;; Should load simple_ec.lisp before this

(defpackage :simple-grammar-based-GP
  (:use :simple-evolutionary-computation :common-lisp)
  (:nicknames :GBGP)
  (:export
   :non-terminal :non-terminal-p :non-terminal-name :non-terminal-rules :non-terminal-min-height
   :rule-p :rule-head :rule-body :rule-min-height :rule-action
   :grammar-non-terminals :grammar-to-internal :pick-rule-for-non-terminal
   :tree-node-p :tree-node-rule :tree-node-children
   :limit-expand-non-terminal-from-grammar :limit-expand-non-terminal
   :for-each-tree-node :for-each-tree-node-depth :parse-tree-nodes :parse-tree-nodes-depth :parse-tree-nodes-height
   :count-parse-tree-nodes :parse-tree-nodes-for-non-terminal
   :parse-tree-sub-node :collect-list :collect-leaf :genotype-to-phenotype
   :*max-parse-tree-height* :*crossover-try-max* :crossover-2-parse-trees :uniform-filter-pick :crossover-2-parse-trees-depth
   :*mutation-max-height* :mutate-parse-tree :mutate-parse-tree-depth
   :steady-state-elitism-GBGP :steady-state-elitism-plain-GBGP))
(in-package :simple-grammar-based-GP)
;;;
;; but for ore efficient lookup of the grammar rules of a non-terminals, the grammar is converted to
;; an internal form.
(defstruct non-terminal
  "One non-terminal.
name: a symbol giving the name of the non-terminal.
rules: a vector of rule structures, giving the alternatives, sorted in ascending order of minimum height.
min-height: minimum parse tree height that can be produced by this non-terminal. 0 for unknown yet."
  name
  rules
  min-height)

(defmethod print-object ((object non-terminal) stream)
  (print-unreadable-object (object stream :type t)
    (format stream "~a :n-rules ~d :min-height ~d"
	    (non-terminal-name object)
	    (length (non-terminal-rules object))
	    (non-terminal-min-height object))))

(defstruct rule
  "One alternative for a non-terminal.
head: head non-terminal for this rule.
body: list of terminals or non-terminals. Anything not a non-terminal is considered a terminal.
min-height: minimum parse tree height that can be produced by this alternative. 0 for unknown yet.
action: binary function of the rule object, and the list of children objects for the non-terminals, should construct the phenotype for this rule."
  head
  body
  min-height
  action)

(defmethod print-object ((object rule) stream)
  (print-unreadable-object (object stream :type t)
    (format stream ":min-height ~d ~a -> " (rule-min-height object) (non-terminal-name (rule-head object)))
    (dolist (x (rule-body object))
      (format stream "~a " (if (non-terminal-p x) (non-terminal-name x) x)))))

(defun grammar-heads (g)
  "Gives the list of non-terminal symbols of the external grammar g."
  (remove-duplicates (mapcar #'car g)))
;;; Common rule-action
(defun collect-list (r cs)
  "non-terminals (tree-node) have been recursively converted, so simply make a copy of the list."
  (declare (ignore r))
  (copy-list cs))
(defun collect-leaf (r cs)
  "cs should contain one item only, take the first."
  (declare (ignore r))
  (car cs))
;;;
(defun rule-alternatives (r)
  (let ((body (cddr r))
	(res nil))
    (loop
	 (let ((next (position :or body)))
	   (unless next
	     (push body res)
	     (return (nreverse res)))
	   (push (subseq body 0 next) res)
	   (setf body (nthcdr (1+ next) body))))))
(defun expand-grammar-short-hands (g)
  "Expand :or and :> for external grammar g, and returns the expanded form of the external grammar."
  (let ((res nil))
    (dolist (r g (nreverse res))
      (let ((nt (first r))
	    (default-action (if (eql (second r) :>) #'collect-leaf #'collect-list)))
	(dolist (a (rule-alternatives r))
	  (push `(,nt -> ,@(if (find :action a)
			       a
			       (append a (list :action default-action))))
		res))))))
;;;
(defun grammar-non-terminals (G)
  "Gives the list of non-terminals of internal grammar G."
  (let ((r nil))
    (maphash #'(lambda (k v) (declare (ignore k)) (push v r)) G)
    r))
(defun calc-grammar-min-heights! (G)
  "Calculate the minimum heights of each non-terminals and each alternative of each non-terminals of the internal grammar G.
Modify G to fill in the min-height fields in the process.
Returns the possibly modified grammar G."
  ; use a naive algorithm to iterative go through all alternatives to determine the minimum heights, until there are no further updates.
  (let ((updatedp nil)
	(ns (grammar-non-terminals G)))
    ; temporarily for the min-height:
    ; 0 means unknown yet.
    ; >0 means min-height known and fully resolved.
    ; note that a min-height of h can be determined in h iterations
    (labels ((calc-rule-min-height! (r)
	       ; update the min-height of rule r, and return the min-height, if the min-heights of
	       ; all its non-terminals in the body is known. Otherwise return 0.
	       (if (zerop (rule-min-height r))
		   (let ((mh 0)
			 (fail? nil)
			 (has-non-terminal? nil))
		     (do ((L (rule-body r) (cdr L)))
			 ((or (null L) fail?)
			  (cond (fail? 0)
				(has-non-terminal?
				 (setf updatedp t)
				 (setf (rule-min-height r) (1+ mh)))
				(t
				 (setf updatedp t)
				 (setf (rule-min-height r) 1))))
		       ;
		       (when (non-terminal-p (car L))
			 (setf has-non-terminal? t)
			 (let ((h (non-terminal-min-height (car L))))
			   (if (zerop h)
			       (setf fail? t)
			       (setf mh (max mh h)))))))
		   (rule-min-height r)))
	     (calc-non-terminal-min-height! (n)
	       ; update the min-height of non-terminal r, and return the min-height if known.
	       ; If the min-height is known, but some alternatives still not resolved, then
	       ; put the min-height to be negative of the min-height.
	       (let ((mh 0))
		 (dolist (r (non-terminal-rules n))
		   (let ((x (calc-rule-min-height! r)))
		     (if (plusp x)
			 (setf mh (if (zerop mh) x (min x mh))))))
		 (when (/= mh (non-terminal-min-height n)) ; different
		   (setf updatedp t)
		   (setf (non-terminal-min-height n) mh)))))
      ;; now iterate, until no changes can be made
      (loop
	 (setf updatedp nil)
	 (dolist (n ns)
	   (calc-non-terminal-min-height! n))
	 (when (not updatedp) (return)))
      ;; done
      G)))

(defun grammar-to-internal (ext-g)
  "Converts grammar into internal form: a hash table (symbol of non-terminals as keys) of non-terminals."
  ;; a hash table using the non-terminals as keys, and the value is the vector of grammar rules for that non-terminal
  (let* ((g (expand-grammar-short-hands ext-g))
	 (tab (make-hash-table :test #'eq))
	 (heads (grammar-heads g)))
    (labels ((new-non-terminal (s)
	       (or (gethash s tab)
		   (setf (gethash s tab) (make-non-terminal :name s :rules nil :min-height 0))))
	     (new-term (s)
	       (if (member s heads)
		   (new-non-terminal s)
		   s))
	     (new-rule (n L act)
	       ; L is list of terminal or non-terminal, for only the rule
	       ; act is the action
	       (make-rule :head n
			  :body (mapcar #'new-term L)
			  :min-height 0
			  :action act)))
      (dolist (r g)
	(let* ((n (new-non-terminal (car r)))
	       (pos (position :action (cddr r)))
	       (rule-ns (if pos
			    (subseq (cddr r) 0 pos)
			    (cddr r)))
	       (act (if pos (nth (1+ pos) (cddr r)) nil)))
	  (push (new-rule n rule-ns act)
		(non-terminal-rules n))))
      ; now get the min-heights right, and sort and convert the rules
      (maphash #'(lambda (k v)
		   (declare (ignore k))
		   (setf (non-terminal-rules v)
			 (coerce (sort (non-terminal-rules v) #'< :key #'rule-min-height) 'vector)))
	       (calc-grammar-min-heights! tab))
      tab)))

;;
(defun pick-rule-for-non-terminal (G s)
  "Gives a random grammar rule from the internal grammar G for non-terminal s (as a symbol)."
  (let ((L (non-terminal-rules (gethash s G))))
    (svref L (random (length L)))))

;;;
;;; An internal node of the parse tree corresponds to a grammar rule of a non-terminal of the grammar,
;;; and anything other than a tree-node is considered a leave in the parse tree.
(defstruct tree-node
  "An internal node of the parse tree.
rule contains the grammar rule corresponding to this node.
children is a list of tree-node for the non-terminals, and the terminals of the grammar rule, in the order they appear in the rule."
  rule
  children)

(defmethod print-object ((object tree-node) stream)
  (print-unreadable-object (object stream :type t)
    (format stream ":rule ~a :children ~a" (tree-node-rule object) (tree-node-children object))))

(defun expand-non-terminals (L expand-f)
  "L is a list of terminals (anything not non-terminal) and non-terminals (the structure). Expand only those non-terminals with expand-f, and leave the terminals unchanged, and return the list of expanded results."
  (cond ((null L) nil)
	((non-terminal-p (car L)) (cons (funcall expand-f (car L)) (expand-non-terminals (cdr L) expand-f)))
	(t (cons (car L) (expand-non-terminals (cdr L) expand-f)))))
(defun random-expand-non-terminal (n)
  "Expand the non-terminal n, until it becomes a parse tree of terminals.
This is dangerous, as no limit to the depth of tree of enforced, so this may cause stackoverflow."
  (let* ((rs (non-terminal-rules n))
	 (r (svref rs (random (length rs)))))
    (make-tree-node
     :rule r
     :children (expand-non-terminals (rule-body r) #'random-expand-non-terminal))))
(defun random-expand-non-terminal-from-grammar (G s)
  "Randomly expand the non-terminal s (as symbol) for the internal grammar G.
This is dangerous, as no limit to the depth of tree of enforced, so this may cause stackoverflow."
  (random-expand-non-terminal (gethash s G)))
;; random expand with limitation in parse tree height
(defun random-limit-rule (rs h)
  "rs is a vector of rules sorted in ascending order of min-height, randomly pick one from those with min-height <= h."
  (do ((L (1- (length rs)) (1- L)))
      ((or (<= L 0) (<= (rule-min-height (svref rs L)) h))
       (svref rs (random (1+ L))))))
(defun limit-expand-non-terminal (n h)
  "Expand the non-terminal n randomly, until it becomes a parse tree of terminals, but limit the total tree height to be <= h."
  (let ((r (random-limit-rule (non-terminal-rules n) h)))
    (make-tree-node
     :rule r
     :children (expand-non-terminals (rule-body r) #'(lambda (z) (limit-expand-non-terminal z (1- h)))))))
(defun limit-expand-non-terminal-from-grammar (G s h)
  "Randomly expand the non-terminal s (as symbol) for the internal grammar G, limiting the initial tree height to be <= h."
  (limit-expand-non-terminal (gethash s G) h))
;;;
;; now the simple crossover and mutation operators.
(defun for-each-tree-node (f tree)
  "Call function f on each of the nodes of tree-node tree and its descendants, in no guaranteed order."
  (funcall f tree)
  (dolist (c (tree-node-children tree))
    (if (tree-node-p c)
	(for-each-tree-node f c))))
(defun for-each-tree-node-depth (f tree)
  "Call the binary function f on each of the nodes of tree-node tree and descendants and its depth, in no guaranteed order.
The root is considered to have depth 0."
  (labels ((traverse (n d)
	     (funcall f n d)
	     (dolist (c (tree-node-children n))
	       (if (tree-node-p c)
		   (traverse c (1+ d))))))
    (traverse tree 0)))
(defun parse-tree-nodes (tree)
  "Gives the list of tree-nodes in the parse-tree tree.
Note that tree-nodes correpsond to non-terminals."
  (let ((ns '()))
    (for-each-tree-node #'(lambda (z) (push z ns)) tree)
    ns))
(defun min-height-of-tree-node (n)
  (non-terminal-min-height (rule-head (tree-node-rule n))))
(defun parse-tree-nodes-depth (tree max-height)
  "Gives the list of pair of tree-nodes that together with its min-height does not exceed the max-height, and its depth in the parse-tree tree.
Note that tree-nodes correspond to non-terminals."
  ;; also considered writing to an adjustable vector (which can be pre-allocated), but keep to this simpler method in fear of premature optimization.
  (let ((ns '()))
    (for-each-tree-node-depth #'(lambda (z d)
				  (if (<= (+ d (min-height-of-tree-node z)) max-height)
				      (push (cons z d) ns)))
			      tree)
    ns))
(defun parse-tree-nodes-height (tree)
  "Gives the list of pair of tree-nodes with its height. Leaf is considered to have height 1."
  (let ((ns nil))
    (labels ((traverse (n)
	       ;; traverse the children of n, and itself, and return its height
	       (let ((h 1))
		 (dolist (c (tree-node-children n))
		   (if (tree-node-p c)
		       (setf h (max h (1+ (traverse c))))))
		 (push (cons n h) ns)
		 h)))
      (traverse tree)
      ns)))
(defun count-parse-tree-nodes (tree)
  (let ((count 0))
    (for-each-tree-node #'(lambda (z) (declare (ignore z)) (incf count)) tree)
    count))
(defun parse-tree-nodes-for-non-terminal (tree n)
  "Gives the list of tree-nodes corresponding to non-terminal n in the parse-tree tree."
  (let ((ns '()))
    (for-each-tree-node #'(lambda (z)
			    (if (eq n (rule-head (tree-node-rule z)))
				(push z ns)))
			tree)
    ns))

(defun parse-tree-sub-node (tree n1 n2)
  "Return a copy of tree, where node n1 is replaced with node n2."
  (if (eq tree n1)
      n2
      (make-tree-node
       :rule (tree-node-rule tree)
       :children (mapcar #'(lambda (z)
			     (if (tree-node-p z)
				 (parse-tree-sub-node z n1 n2)
				 z))
			 (tree-node-children tree)))))
;;;
(defun genotype-to-phenotype (tree)
  "Use the actions of each tree-node (if non-nil) to recursively construct the phenotype of the parse tree.
If no action is given for the rule, the tree-node itself is returned."
  (let ((act (or (rule-action (tree-node-rule tree))
		 #'collect-list)))
    (funcall act
	     (tree-node-rule tree)
	     (mapcar #'(lambda (c) (if (tree-node-p c)
				       (genotype-to-phenotype c)
				       c))
		     (tree-node-children tree)))))
;;
(defparameter *max-parse-tree-height* 10)
;;
(defparameter *crossover-try-max* 10)
(defun crossover-2-parse-trees (t1 t2)
  "Crossover two parse trees (tree-node) t1 and t2 by randomly picking a tree-node (for non-terminal) in t1, and randomly picking a tree-node for the same non-terminal in t2, and return a copy of t1 with the chosen node replaced with that from t2.
Note that t1 and t2 should be generated from the same instance of internal grammar, so that the same non-terminal are the same object instance.
If after *crossover-try-max* times of trying, and still cannot randomly pick a non-terminal in t1 that occurs in t2, then returns t1 unchanged."
  (let* ((ns (parse-tree-nodes t1))
	 (L (length ns)))
    (do ((i *crossover-try-max* (1- i)))
	((<= i 0) t1)
      (let* ((n1 (nth (random L) ns))
	     (nt1 (rule-head (tree-node-rule n1)))
	     (ns2 (parse-tree-nodes-for-non-terminal t2 nt1))
	     (n2 (if (null ns2) nil (nth (random (length ns2)) ns2))))
	(if n2 (return (parse-tree-sub-node t1 n1 n2)))))))
(defun uniform-filter-pick (lst test)
  "Among those items in lst that satify the unary function test, uniformly pick one. Use an online method, so need not create a temporary list, but will use more random numbers. Gives nil if no item satify the condition."
  (let ((r nil)
	(k 1))
    (dolist (c lst r)
      (when (funcall test c)
	(if (or (= k 1) (= (random k) 1))
	    (setf r c))
	(incf k)))))
(defun crossover-2-parse-trees-depth (t1 t2)
  "Crossover two parse trees (tree-node) t1 and t2 by randomly picking a tree-node (for non-terminal) in t1, and randomly picking a tree-node for the same non-terminal in t2, and return a copy of t1 with the chosen node replaced with that from t2.
Note that t1 and t2 should be generated from the same instance of internal grammar, so that the same non-terminal are the same object instance.
If after *crossover-try-max* times of trying, and still cannot randomly pick a non-terminal in t1 that occurs in t2, then returns t1 unchanged.
Will try to maintain that the tree depth does not exceed *max-parse-tree-height*"
  (let* ((ns1 (parse-tree-nodes-depth t1 *max-parse-tree-height*)) ; ((node . depth) ...)
	 (ns2 (parse-tree-nodes-height t2)) ; ((node . height) ...)
	 (L (length ns1)))
    (if (null ns1)
	t1
	(do ((i *crossover-try-max* (1- i)))
	    ((<= i 0) t1)
	  (let* ((nr (nth (random L) ns1)) ; nr is (node . depth)
		 (n1 (car nr))
		 (d1 (cdr nr))
		 (nt1 (rule-head (tree-node-rule n1)))
		 (n2 (car (uniform-filter-pick ns2 #'(lambda (z) (and (eq (rule-head (tree-node-rule (car z))) nt1)
								      (<= (+ d1 (cdr z)) *max-parse-tree-height*)))))))
	    (if n2 (return (parse-tree-sub-node t1 n1 n2))))))))
;;
; the maximum height of random subtree, but will use the min-height of the non-terminal being replaced if it is larger.
(defparameter *mutation-max-height* 5)
(defun mutate-parse-tree (tree)
  "Randomly pick a tree-node from the parse tree, and replace it with a random subtree from the same non-terminal."
  (let* ((ns (parse-tree-nodes tree))
	 (n (nth (random (length ns)) ns)) ; random tree-node
	 (nt (rule-head (tree-node-rule n))) ; the non-terminal
	 (new-n (limit-expand-non-terminal nt (max *mutation-max-height* (non-terminal-min-height nt)))))
    (parse-tree-sub-node tree n new-n)))
(defun mutate-parse-tree-depth (tree)
  "Randomly pick a tree-node from the parse tree, and replace it with a random subtree from the same non-terminal.
Will try to maintain that the tree depth does not exceed *max-parse-tree-height*"
  (let* ((ns (parse-tree-nodes-depth tree *max-parse-tree-height*))
	 (nr (if (null ns)
		 (cons tree 0)
		 (nth (random (length ns)) ns))) ; random (tree-node . depth)
	 (n (car nr))
	 (d (cdr nr))
	 (nt (rule-head (tree-node-rule n))) ; the non-terminal
	 (new-n (limit-expand-non-terminal nt (max (min *mutation-max-height*
							(- *max-parse-tree-height* d))
						   (non-terminal-min-height nt)))))
    (parse-tree-sub-node tree n new-n)))
;;
;;;
;; Copied and modified from generational-elitism-EC to use GBGP
;; This can be a simple template for further customization.
(defun steady-state-elitism-GBGP (internal-grammar start-symbol evaluator
				  &key (population-size 500) (p-crossover 0.9) (p-mutation 0.05)
				  (generations 50) (init-max-height 5) (max-tree-height 10)
				    (max-evals nil) (optimum nil) (stop-at-optimum nil)
				  (better-score #'>) (report-fitness-record nil))
  ;; evaluator should give a fitness of a chromosome
  ;; (better-score a b) says score a given by evaluator is better than score b
  ;; also record the fitness changes
  ;; uses 2-tournament
  ;; returns the best found (phenotype), best found (genotype), number of evaluations used, and the last generation
  (let ((n-evals 0))
    (labels ((random-parse-tree ()
	       (limit-expand-non-terminal-from-grammar internal-grammar start-symbol init-max-height))
	     (tree-evaluator (chr)
	       (incf n-evals)
	       (cons (funcall evaluator chr)
		     (count-parse-tree-nodes chr)))
	     (better-fitness (ind1 ind2)
	       ;; individual ind1 better than ind2
	       ;; here the fitnesses are (score . tree-size)
	       (let ((f1 (individual-fitness ind1))
		     (f2 (individual-fitness ind2)))
		 (or (funcall better-score (car f1) (car f2))
		     (and (not (funcall better-score (car f2) (car f1)))
			  (< (cdr f1) (cdr f2))))))
	     (better-individual (ind1 ind2)
	       (if (better-fitness ind1 ind2) ind1 ind2))
	     (ind-crossoveror (ind1 ind2)
	       (new-individual (crossover-2-parse-trees-depth (individual-chromosome ind1)
							      (individual-chromosome ind2))
			       #'tree-evaluator))
	     (ind-mutator (ind)
	       (new-individual (mutate-parse-tree-depth (individual-chromosome ind))
			       #'tree-evaluator)))
      (let* ((*max-parse-tree-height* max-tree-height)
	     (*mutation-max-height* max-tree-height)
	     (n-crossover (ceiling (* population-size p-crossover)))
	     (n-mutation (ceiling (* population-size p-mutation)))
	     (n-reproduction (max 0 (- population-size n-crossover n-mutation 1)))
	     (cur-pop (fill-all-with (make-array population-size)
			(new-individual (random-parse-tree) #'tree-evaluator)))
	     (best-so-far (best-of cur-pop #'better-fitness))
	     (next-pop (make-array population-size :initial-element nil))
	     (selector (new-2tournament-selector #'better-individual))
	     (crossoveror (if report-fitness-record
			      (add-fitness-record-crossover #'ind-crossoveror)
			      #'ind-crossoveror))
	     (mutator (if report-fitness-record
			  (add-fitness-record-mutation #'ind-mutator)
			  #'ind-mutator)))
	;;
	(dotimes (g generations)
	  ;; next generation
	  (let ((sel (tournament-prepare cur-pop))
		(cur-best (best-of cur-pop #'better-fitness)))
	    (fill-with next-pop
	      (1 cur-best)
	      (n-crossover (funcall crossoveror (funcall selector sel) (funcall selector sel)))
	      (n-mutation (funcall mutator (funcall selector sel)))
	      (n-reproduction (funcall selector sel)))
	    ;; sort by fitness
	    (sort next-pop #'better-fitness)
	    ;; report
	    (if (better-fitness cur-best best-so-far) (setf best-so-far cur-best))
	    (format t "===Generation ~d ===~%" (1+ g))
	    (format t "~%best-so-far: fitness: ~a " (individual-fitness best-so-far))
	    (princ (genotype-to-phenotype (individual-chromosome best-so-far)))
	    (terpri)
	    ;; fitness record
	    (when report-fitness-record
	      (report-fitness-changes (1+ g))))
	  ;; terminate early
	  (when (or (and max-evals (>= n-evals max-evals))
		    (and stop-at-optimum
			 (equal optimum (car (individual-fitness best-so-far)))))
	    (return))
	  ;;
	  (let ((tmp cur-pop))
	    (setf cur-pop next-pop)
	    (setf next-pop tmp)))
	;;
	(values (genotype-to-phenotype (individual-chromosome best-so-far)) best-so-far n-evals cur-pop)))))

;;;
;;; without using tree size in selection
(defun steady-state-elitism-plain-GBGP (internal-grammar start-symbol evaluator
					&key (population-size 500) (p-crossover 0.9) (p-mutation 0.05)
					  (generations 50) (init-max-height 5) (max-tree-height 10)
					  (max-evals nil) (optimum nil) (stop-at-optimum nil)
					  (better-score #'>) (report-fitness-record nil))
  ;; evaluator should give a fitness of a chromosome
  ;; (better-score a b) says score a given by evaluator is better than score b
  ;; also record the fitness changes
  ;; uses 2-tournament
  ;; returns the best found (phenotype), best found (genotype), number of evaluations used, and the last generation
  (let ((n-evals 0))
    (labels ((random-parse-tree ()
	       (limit-expand-non-terminal-from-grammar internal-grammar start-symbol init-max-height))
	     (tree-evaluator (chr)
	       (incf n-evals)
	       (funcall evaluator chr))
	     (better-fitness (ind1 ind2)
	       ;; individual ind1 better than ind2
	       ;; here the fitnesses are (score . tree-size)
	       (let ((f1 (individual-fitness ind1))
		     (f2 (individual-fitness ind2)))
		 (funcall better-score f1 f2)))
	     (better-individual (ind1 ind2)
	       (if (better-fitness ind1 ind2) ind1 ind2))
	     (ind-crossoveror (ind1 ind2)
	       (new-individual (crossover-2-parse-trees-depth (individual-chromosome ind1)
							      (individual-chromosome ind2))
			       #'tree-evaluator))
	     (ind-mutator (ind)
	       (new-individual (mutate-parse-tree-depth (individual-chromosome ind))
			       #'tree-evaluator)))
      (let* ((*max-parse-tree-height* max-tree-height)
	     (*mutation-max-height* max-tree-height)
	     (n-crossover (ceiling (* population-size p-crossover)))
	     (n-mutation (ceiling (* population-size p-mutation)))
	     (n-reproduction (max 0 (- population-size n-crossover n-mutation 1)))
	     (cur-pop (fill-all-with (make-array population-size)
			(new-individual (random-parse-tree) #'tree-evaluator)))
	     (best-so-far (best-of cur-pop #'better-fitness))
	     (next-pop (make-array population-size :initial-element nil))
	     (selector (new-2tournament-selector #'better-individual))
	     (crossoveror (if report-fitness-record
			      (add-fitness-record-crossover #'ind-crossoveror)
			      #'ind-crossoveror))
	     (mutator (if report-fitness-record
			  (add-fitness-record-mutation #'ind-mutator)
			  #'ind-mutator)))
	;;
	(dotimes (g generations)
	  ;; next generation
	  (let ((sel (tournament-prepare cur-pop))
		(cur-best (best-of cur-pop #'better-fitness)))
	    (fill-with next-pop
	      (1 cur-best)
	      (n-crossover (funcall crossoveror (funcall selector sel) (funcall selector sel)))
	      (n-mutation (funcall mutator (funcall selector sel)))
	      (n-reproduction (funcall selector sel)))
	    ;; sort by fitness
	    (sort next-pop #'better-fitness)
	    ;; report
	    (if (better-fitness cur-best best-so-far) (setf best-so-far cur-best))
	    (format t "===Generation ~d ===~%" (1+ g))
	    (format t "~%best-so-far: fitness: ~a " (individual-fitness best-so-far))
	    (princ (genotype-to-phenotype (individual-chromosome best-so-far)))
	    (terpri)
	    ;; fitness record
	    (when report-fitness-record
	      (report-fitness-changes (1+ g))))
	  ;; terminate early
	  (when (or (and max-evals (>= n-evals max-evals))
		    (and stop-at-optimum
			 (equal optimum (individual-fitness best-so-far))))
	    (return))
	  ;;
	  (let ((tmp cur-pop))
	    (setf cur-pop next-pop)
	    (setf next-pop tmp)))
	;;
	(values (genotype-to-phenotype (individual-chromosome best-so-far)) best-so-far n-evals cur-pop)))))

;;;
