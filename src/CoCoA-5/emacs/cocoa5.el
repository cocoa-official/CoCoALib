;; cocoa5.el --- major mode for editing and running CoCoA

;; Copyright (c) 2005-2019  John Abbott, Anna Bigatti
;; This file is part of the CoCoA distribution.
;; You are free to use any part of this file in your own programs.
;; Many thanks to Giovanna Roda and Giorgio Vecchiocattivi 
;;   for their substantial early contributions.
;; cocoa-comint-mode introduced by Burkhard Zimmermann (May 2005)

;; Jun 2005: 
;; tab does both indentation and dabbrev completion -- B.Zimmermann
;; C-c C-e for EndDefine also prints comment with function name
;; F2 goes to error line after sourcing a file -- B.Zimmermann

;; Oct 2005: 
;; "CoCoA not responding?" (panic button) in the menu
;; "Source File Parse Error Line" in the menu

;; Nov 2005: 
;; toggle "TAB for Completion" in the menu

;; May 2009
;; Fixed a bug which prevented recognition of Define at the start of a file.
;; Minor code cleaning.

;; July 2010
;; Fixed a bug on escaped quotes inside strings

;; May 2011
;; Fixed abbrev in comments and strings!
;; fixed call to abbrev when typing semicolon
;; removed semicolon fro abbrev expansions

;; June 2012
;; added customization for --fullCoCoALibError

;; November 2012
;; added C-c ]  for End..
;; fixed wrong indentation within parentheses

;; March 2013
;; new code for send-region (calling SourceRegion)
;; and disabled echoing of the region during SourceRegion.

;; June 2013 (EACA school Valladolid)
;; added C-RET for sending line
;; changed C-c C-m into C-c C-h for calling manual
;; some improvements to the menu

;; April 2014 (for Easter release)
;; recognizing lower case keywords
;; menu for toggling automatic capitalization

;; May 2014
;; syntax colours also in *cocoa5* buffer
;; prompt can be configured via emacs
;; path corrections for windows ("\" no longer needed)

;; August 2014
;; Added cocoa5-send-line-or-region; bindings to C-RET M-RET
;; menu: side-by-side/top-bottom

;; April 2017
;; send interrupt to CoCoA (C-c C-c) -- also in menu

;; November 2019
;; Added menu option to allow cocoa5 to be started with --enable-system-cmd

;; May 2020
;; long lines & lines containing "strange" chars will be sent via send-region

;; March 2022
;; help for word also for comint
;; cocoa5-close-block also for comint

(setq cocoa5-mode-version "07 March 2022")

;;==============================================================>
;;     CUSTOMIZATION in cocoa5.emacs
;;==============================================================>

;;==============================================================>
;  \M-x cocoa5
;; will open a buffer called *cocoa5* with a running CoCoA

;; you can edit a CoCoA file and send commands from it to *cocoa5*:
;   "\C-c\C-l" 'cocoa5-send-line
;   "\C-c\C-r" 'cocoa5-send-region
;   "\C-c\C-f" 'cocoa5-source-file
;   "\C-c\C-c" 'cocoa5-interrupt-or-comment-region
;   "\C-c\C-e" 'cocoa5-close-block  ;; writes the appropriate "End"
;   "\C-c ]"   'cocoa5-close-block  ;; same as above
;   "\C-c\C-w" 'cocoa5-open-wordlist-file
;   "\C-c\C-q" 'cocoa5-shell-quit
;   "\C-c\C-h" 'cocoa5-word-help
;;==============================================================>


(defgroup cocoa5 nil
  "Options for cocoa mode.  'Font Lock' to change syntax colours"
  :prefix "cocoa5-"
  :link '(custom-group-link :tag "Font Lock Faces group" font-lock-faces)
  :group 'languages)

(defcustom cocoa5-tab-is-completion nil
  "Non-nil means that TAB runs auto-completion\n\
   This is in addition to  M-/  which runs auto-completion anyway.\n\
   Effect is immediate."
  :type 'boolean
  :group 'cocoa5)

(defcustom cocoa5-split-is-side-by-side nil
  "Non-nil means that source/shell are side/side instead of top/bottom\n\
   Effect is immediate."
  :type 'boolean
  :group 'cocoa5)

(defcustom cocoa5-enable-system-cmd-flag nil
  "Non-nil means that a new cocoa5 session will accept SystemCommand\n\
   Must be set BEFORE starting cocoa5."
  :type 'boolean
  :group 'cocoa5)

(defcustom cocoa5-echo-source-region nil
  "CoCoA-5 SourceRegion verboseness:\n\
   Non-nil means that SourceRegion echoes the read region as a comment.\n\
   Effect is immediate."
  :type 'boolean
  :group 'cocoa5)

(defcustom cocoa5-warning-level "-w2"
  "Warning level for CoCoA-5:\n\
   -w0 none, -w1 low, -w2 normal (the default), -w3 pedantic\n\
   of course, the higher the setting, the better the code!\n\
   Effect when starting a new CoCoA-5."

  :type 'string
  :group 'cocoa5)

(defcustom cocoa5-prompt cocoa5-default-prompt
  "Prompt for CoCoA-5.\n\
   Effect when starting a new CoCoA-5"
  :type 'string
  :group 'cocoa5)

(defcustom cocoa5-full-cocoalib-error  nil
  "CoCoA-5 errors verboseness:\n\
   Enables a verbose error reporting from CoCoALib,\n\
   useful if you want to know more about the errors you get.\n\
   Effect when starting a new CoCoA-5"
  :type 'boolean
  :group 'cocoa5)

;;; End customizable variables.


;;-------------------------------------------------------------;;
;;                    EDITING CoCoA files                      ;;
;;-------------------------------------------------------------;;

;;-------------------------------------------------------------->
;; cocoa5-mode-abbrev-table
;;-------------------------------------------------------------->
;; We use it to capitalize the CoCoA keywords .

(defvar cocoa5-mode-abbrev-table nil
  "Abbrev table in use in CoCoA buffers.")

(define-abbrev-table 'cocoa5-mode-abbrev-table '(
    ("alias"   		"Alias"   	cocoa5-check-expansion 0)
    ("and"    		"And"   	cocoa5-check-expansion 0)
    ("break"   		"Break"   	cocoa5-check-expansion 0)
    ("ciao"    		"Ciao"   	cocoa5-check-expansion 0)
    ("continue"   	"Continue"   	cocoa5-check-expansion 0)
    ("define"  		"Define"  	cocoa5-check-expansion 0)
    ("describe"        	"Describe"  	cocoa5-check-expansion 0) ;; JAA: added 2017-12-15
    ("deglex"		"DegLex"     	cocoa5-check-expansion 0)
    ("degrevlex"	"DegRevLex"     cocoa5-check-expansion 0)
    ("do"      		"Do"      	cocoa5-check-expansion 0)
    ("elif"    		"Elif"    	cocoa5-check-expansion 0)
    ("else"    		"Else"    	cocoa5-check-expansion 0)
    ("end"     		"End"    	cocoa5-check-expansion 0)
    ("endblock"         "EndBlock"   	cocoa5-check-expansion 0)
    ("endcatch"         "EndCatch"   	cocoa5-check-expansion 0)
    ("enddefine"        "EndDefine"   	cocoa5-check-expansion 0)
    ("endfor"           "EndFor"   	cocoa5-check-expansion 0)
    ("endforeach"       "EndForeach"  	cocoa5-check-expansion 0)
    ("endif"            "EndIf"   	cocoa5-check-expansion 0)
    ("endfunc"          "EndFunc"   	cocoa5-check-expansion 0)
    ("endpackage"       "EndPackage"  	cocoa5-check-expansion 0)
    ("endusing"         "EndUsing"   	cocoa5-check-expansion 0)
    ("endwhile"         "EndWhile"   	cocoa5-check-expansion 0)
    ("exit"    		"Exit"   	cocoa5-check-expansion 0)
    ("false"     	"False"     	cocoa5-check-expansion 0)
    ("for"     		"For"     	cocoa5-check-expansion 0)
    ("foreach" 		"Foreach" 	cocoa5-check-expansion 0)
    ("func" 		"Func"          cocoa5-check-expansion 0)
    ("ideal"      	"ideal"      	cocoa5-check-expansion 0)
    ("if"      		"If"      	cocoa5-check-expansion 0)
    ("in"      		"In"      	cocoa5-check-expansion 0)
    ("isin"            	"IsIn"      	cocoa5-check-expansion 0)
;;    ("indentation"	"Indentation"	cocoa5-check-expansion 0) ;; JAA: removed 2014-03-11
    ("lex"    		"Lex"  	        cocoa5-check-expansion 0)
    ("mat"    		"mat"  	        cocoa5-check-expansion 0)
    ("module"    	"Module"        cocoa5-check-expansion 0)
    ("on"    		"On"   	        cocoa5-check-expansion 0)
    ("or"    		"Or"   	        cocoa5-check-expansion 0)
    ("package"          "Package"  	cocoa5-check-expansion 0)
    ("posto"   		"PosTo"   	cocoa5-check-expansion 0)
    ("print"   		"Print"   	cocoa5-check-expansion 0)
    ("println" 		"PrintLn" 	cocoa5-check-expansion 0)
    ("protect" 		"Protect" 	cocoa5-check-expansion 0) ;; JAA: added 2017-12-15
    ("quit"    		"Quit"   	cocoa5-check-expansion 0)
    ("record"  		"Record"  	cocoa5-check-expansion 0)
    ("return"  		"Return"  	cocoa5-check-expansion 0)
;;    ("set"     		"Set"   	cocoa5-check-expansion 0) ;; JAA: removed 2014-03-11
    ("source"    	"Source"   	cocoa5-check-expansion 0) ;; JAA: added 2017-12-15
    ("sourceregion"    	"SourceRegion" 	cocoa5-check-expansion 0) ;; JAA: added 2017-12-15
    ("step"    		"Step"   	cocoa5-check-expansion 0)
    ("toplevel"    	"TopLevel"   	cocoa5-check-expansion 0) ;; AMB: added 2010-07-19
    ("quit"    		"Quit"   	cocoa5-check-expansion 0)
    ("then"    		"Then"    	cocoa5-check-expansion 0)
;;    ("time"    		"Time"    	cocoa5-check-expansion 0) ;; JAA: removed 2014-03-11
    ("true"    		"True"    	cocoa5-check-expansion 0)
    ("to"      		"To"      	cocoa5-check-expansion 0)
    ("topos"   		"ToPos"      	cocoa5-check-expansion 0)
    ("try"   		"Try"      	cocoa5-check-expansion 0) ;; AMB: added 2010-08-03
;;    ("unset"     	"Unset"     	cocoa5-check-expansion 0) ;; JAA: removed 2014-03-11
    ("unprotect"     	"Unprotect"    	cocoa5-check-expansion 0) ;; JAA: added 2017-12-15
    ("until"     	"Until"     	cocoa5-check-expansion 0) ;; JAA: added 2012-04-30
    ("uponerror"	"UponError"    	cocoa5-check-expansion 0) ;; AMB: added 2010-08-03
    ("use"     		"Use"     	cocoa5-check-expansion 0)
;;    ("using"   		"Using"   	cocoa5-check-expansion 0) ;; JAA: removed 2014-03-11
;;    ("var"   		"Var"   	cocoa5-check-expansion 0) ;; JAA: removed 2014-03-11
    ("vector"  		"Vector"   	cocoa5-check-expansion 0)
    ("while"   		"While"   	cocoa5-check-expansion 0)
    ("xel"    		"Xel"  	        cocoa5-check-expansion 0)
    ))


(defun cocoa5-inside-comment-or-string-p ()
  "Check if the point is inside a comment or a string."
  (save-excursion
    (let* ((origpoint (point))
           state)
      (goto-char 1)
      (while (> origpoint (point))
	(setq state (parse-partial-sexp (point) origpoint 0)))
      (or (nth 3 state) (nth 4 state))))  )  ;; 4 comment

(defun cocoa5-check-expansion ()
  "If abbrev was made within a comment or a string, de-abbrev!"
  (if (cocoa5-inside-comment-or-string-p)
	(unexpand-abbrev)))

(defvar cocoa5-imenu-generic-expression
  '("^[ \t]*\\([Dd]efine\\)[ \t\n]+\\([a-zA-Z0-9_.:]+\\)" . (2))
  "Imenu expression for CoCoA-mode.  See `imenu-generic-expression'.")

;; JAA removed Alias & Cond  20090512
;; JAA removed Repeat 20091103 && reinserted in 20120127
;; AMB removed Using  20140404

(defconst cocoa5-lc-subblock
  "block\\|catch\\|for\\|foreach\\|func\\|if\\|repeat\\|try\\|while" )
(defconst cocoa5-uc-subblock
  "Block\\|Catch\\|For\\|Foreach\\|Func\\|If\\|Repeat\\|Try\\|While" )

(defconst cocoa5-lc-startkeywords 
  (concat cocoa5-lc-subblock "\\|define\\|package" ))
(defconst cocoa5-uc-startkeywords
  (concat cocoa5-uc-subblock "\\|Define\\|Package" ))
(defconst cocoa5-startkeywords-list
  (concat "\\(" cocoa5-lc-startkeywords "\\|" cocoa5-uc-startkeywords "\\)" ))

(defconst cocoa5-endkeywords-list 
  (concat "\\([Uu]ntil\\)\\|"
	  "\\("  (concat "end\\(" cocoa5-lc-startkeywords "\\)")
	  "\\|"  (concat "End\\(" cocoa5-uc-startkeywords "\\)")
	  "\\)"))

;;-------------------------------------------------------------->
;; Regular expressions used to calculate indent, etc.
;;-------------------------------------------------------------->

(defconst cocoa5-symbol-re      "\\<[a-zA-Z_][a-zA-Z_0-9.]*\\>")
;; (defconst cocoa5-beg-block-re ;; JAA removed Alias & Cond  20090512
(defconst cocoa5-beg-block-re (concat "\\<" cocoa5-startkeywords-list "\\>"))
(defconst cocoa5-end-block-re
  (concat "\\<\\("
	  cocoa5-endkeywords-list
	  "\\)\\>"))
(defconst cocoa5-declaration-re "\\<\\(NIENTE_DA_FISSARE\\)\\>")
(defconst cocoa5-defun-re       "\\<\\([Dd]efine\\)\\>")
(defconst cocoa5-subblock-re
  (concat "\\<\\(" cocoa5-lc-subblock "\\|" cocoa5-uc-subblock "\\)\\>"))

(defconst cocoa5-noindent-re
(concat "\\<\\("
	cocoa5-endkeywords-list "\\|[Ee]xport\\|"
	"[Ee]lse\\|[Ee]lif\\|UponError\\|uponerror" ;; JAA replaced Elsif by Elif 20090512
	"\\)\\>"))
;; ;; (defconst cocoa5-nosemi-re
;;   "\\<\\([Tt]hen\\|[Dd]o\\|[Ee]lse\\|UponError\\|uponerror\\)\\>")
(defconst cocoa5-autoindent-lines-re
  (concat "\\<\\("
	  cocoa5-startkeywords-list
	  "\\|" cocoa5-noindent-re
	  "\\)\\>"))

;; ;;; Strings used to mark beginning and end of excluded text
;; (defconst cocoa5-exclude-str-start "{-----\\/----- EXCLUDED -----\\/-----")
;; (defconst cocoa5-exclude-str-end " -----/\\----- EXCLUDED -----/\\-----}")

;;-------------------------------------------------------------->
;; cocoa5-mode-syntax-table
;;-------------------------------------------------------------->

(defvar cocoa5-mode-syntax-table nil
  "Syntax table in use in cocoa5-mode buffers.")

(if cocoa5-mode-syntax-table
    ()
  (setq cocoa5-mode-syntax-table (make-syntax-table))
  (modify-syntax-entry ?\( "()"  cocoa5-mode-syntax-table)
  (modify-syntax-entry ?\) ")("  cocoa5-mode-syntax-table)
;  (modify-syntax-entry ?* "." cocoa5-mode-syntax-table)
  (modify-syntax-entry ?+ "."    cocoa5-mode-syntax-table)
  (modify-syntax-entry ?= "."    cocoa5-mode-syntax-table)
  (modify-syntax-entry ?% "."    cocoa5-mode-syntax-table)
  (modify-syntax-entry ?< "."    cocoa5-mode-syntax-table)
  (modify-syntax-entry ?> "."    cocoa5-mode-syntax-table)
  (modify-syntax-entry ?& "."    cocoa5-mode-syntax-table)
  (modify-syntax-entry ?| "."    cocoa5-mode-syntax-table)
  (modify-syntax-entry ?_ "w"    cocoa5-mode-syntax-table)
  ; string
;  (modify-syntax-entry ?\' "\""  cocoa5-mode-syntax-table) ;; AMB: rm 20111012
  ; testing paired delimiters
  ; (modify-syntax-entry ?$ "$"    cocoa5-mode-syntax-table)
  ; -- comments 
  (modify-syntax-entry ?- ". 12"    cocoa5-mode-syntax-table)
  ; // and /* */ comments 
  (modify-syntax-entry ?/  ". 124" cocoa5-mode-syntax-table)
  (modify-syntax-entry ?*  ". 23b" cocoa5-mode-syntax-table)
  (modify-syntax-entry ?\n ">"   cocoa5-mode-syntax-table)
)

;;-------------------------------------------------------------->
;; font-lock-keywords
;;-------------------------------------------------------------->

(defvar cocoa5-font-lock-keywords
  (list
   ; testing:  ] first character; - first or after ]
;;    '("\\*\\*\\*[^;]+\\*\\*\\*" . 'font-lock-string-face)
    '("\\*\\*\\*" . 'font-lock-constant-face)
   ; function definition
   '("^[ \t]*\\([Dd]efine\\)\\>[ \t]*\\(\\sw+\\)?"
     (1 font-lock-keyword-face) (2 font-lock-function-name-face nil t))
   ; package definition
   '("^[ \t]*\\([Pp]ackage\\)\\>[ \t]*\\(\\sw+\\)?"
     (1 font-lock-keyword-face) (2 font-lock-function-name-face nil t))
   ; types
   (cons    (concat "\\<\\("
	   "BOOL\\|ERROR\\|FUNCTION\\|"
	   "IDEAL\\|INT\\|LIST\\|"
	   "MAT\\|MATRIXROW\\|MODULE\\|MODULEELEM\\|"
	   "OSTREAM\\|PACKAGE\\|" ;; AMB 2010-08-25 (rm NULL,DEVICE)
	   "RAT\\|RECORD\\|RING\\|RINGELEM\\|RINGHOM\\|"
	   "STRING\\|TAGGED\\|TYPE\\|VOID"
	   "\\)\\>")
	 'font-lock-type-face)
   ; constants
   (cons   (concat "\\<\\("
		   "TRUE\\|FALSE\\|True\\|False\\|true\\|false\\|"
		   "[Ll]ex\\|[Xx]el\\|DegLex\\|DegRevLex\\|"
		   "PosTo\\|ToPos\\)\\>")
	   'font-lock-constant-face)
   ; Ideal, Module...  AMB 2014-07 added "not"
   (cons   (concat "\\<\\("
		   "Ideal\\|Mat\\|Not\\|Record\\|Error\\|"
		   "ideal\\|mat\\|matrix\\|not\\|record\\|error\\|submodule"
		   "\\)\\>")
	   'font-lock-builtin-face)
   ; keywords
   (concat "\\<\\("
	   cocoa5-startkeywords-list "\\|"
	   cocoa5-noindent-re "\\|"
	   "[Aa]nd\\|[Oo]r\\|"
	   "[Bb]reak\\|[Cc]ontinue\\|[Cc]iao\\|[Dd]o\\|[Dd]escribe\\|[Ee]xit\\|"
           "ImportByValue\\|ImportByRef\\|importbyvalue\\|importbyref\\|"
	   "[Ii]n\\|isin\\|IsIn\\|"
	   "[Oo]n\\|[Oo]pt\\|Print\\(\\|Ln\\)\\|print\\(\\|ln\\)\\|[Pp]rotect\\|"
	   "[Qq]uit\\|[Rr]ef\\|[Rr]eturn\\|[Ss]tep\\|" ;; JAA removed Set 20140318
           "Source\\(\\|Region\\)\\|source\\(\\|region\\)\\|" ;; JAA added 2017-12-15
	   "toplevel\\|TopLevel\\|" ;; AMB: added 2010-07-19
	   "[Tt]hen\\|[Tt]o\\|"  ;; JAA removed Time 20140318
	   "[Uu]nprotect\\|[Uu]se" ;; JAA removed Unset 20140318
	   "\\)\\>")
   '("\\<\\(NULLA_DA_FISSARE\\)\\>[ \t]*\\([0-9]+\\)?"
     (1 font-lock-keyword-face) (2 font-lock-preprocessor-face nil t))
   )
  "Additional expressions to highlight in CoCoA mode.")

;;-------------------------------------------------------------->

(defvar cocoa5-indent-level 2
  "*Indentation of CoCoA statements with respect to containing block.")

(defvar cocoa5-case-indent 2
  "*Indentation for case statements.")

(defvar cocoa5-auto-newline nil
  "*Non-nil means automatically newline after semicolons and the punctuation
mark after an end.")

(defvar cocoa5-tab-always-indent t
  "*Non-nil means TAB in CoCoA mode should always reindent the current line,
regardless of where in the line point is when the TAB command is used.")

(defvar cocoa5-auto-endcomments t
  "*Non-nil means a comment { ... } is set after the ends which ends cases and
functions. The name of the function or case will be set between the braces.")

(defvar cocoa5-auto-lineup '(all)
  "*List of contexts where auto lineup of :'s or ='s should be done.
Elements can be of type: 'paramlist', 'declaration' or 'case', which will
do auto lineup in parameterlist, declarations or case-statements
respectively. The word 'all' will do all lineups. '(case paramlist) for
instance will do lineup in case-statements and parameterlist, while '(all)
will do all lineups.")

;; (defvar cocoa5-toggle-completions nil
;;   "*Non-nil means \\<cocoa5-mode-map>\\[cocoa5-complete-word] should try all
;; possible completions one by one.
;; Repeated use of \\[cocoa5-complete-word] will show you all of them.
;; Normally, when there is more than one possible completion,
;; it displays a list of all possible completions.")

;;;
;;;  Macros
;;;

(defsubst cocoa5-get-beg-of-line (&optional arg)
  (save-excursion
    (beginning-of-line arg)
    (point)))

(defsubst cocoa5-get-end-of-line (&optional arg)
  (save-excursion
    (end-of-line arg)
    (point)))

(defun cocoa5-declaration-end ()
  (let ((nest 1))
    (while (and (> nest 0)
		(re-search-forward
		 "[:=]\\|\\(\\<record\\>\\)\\|\\(\\<end\\>\\)"
		 (save-excursion (end-of-line 2) (point)) t))
      (cond ((match-beginning 1) (setq nest (1+ nest)))
	    ((match-beginning 2) (setq nest (1- nest)))
	    ((looking-at "[^(\n]+)") (setq nest 0))))))


(defun cocoa5-declaration-beg ()
  (let ((nest 1))
    (while (and (> nest 0)
		(re-search-backward
"[:=]\\|\\<\\(NULLA_DECLARATION_BEG\\)\\>\\|\\(\\<record\\>\\)\\|\\(\\<end\\>\\)
" (cocoa5-get-beg-of-line 0) t))
      (cond ((match-beginning 1) (setq nest 0))
	    ((match-beginning 2) (setq nest (1- nest)))
	    ((match-beginning 3) (setq nest (1+ nest)))))
    (= nest 0)))


(defsubst cocoa5-within-string ()
  (save-excursion
    (nth 3 (parse-partial-sexp (cocoa5-get-beg-of-line) (point)))))




;;======================================================================
;;;  Electric functions
;;======================================================================
(defun cocoa5-electric-terminate-line ()
  "Terminate line and indent next line."
  (interactive)
  ;; First, check if current line should be indented
  (save-excursion
    (beginning-of-line)
    (skip-chars-forward " \t")
    (if (looking-at cocoa5-autoindent-lines-re)  (cocoa5-indent-line))
    )
  (delete-horizontal-space) ; Removes trailing whitespaces
  (newline)
  ;; Indent next line
  (cocoa5-indent-line)
  ;; Maybe we should set some endcomments
  (if cocoa5-auto-endcomments  (cocoa5-set-auto-comments))
  ;; Check if we shall indent inside comment
  (let ((setstar nil))
    (save-excursion
      (forward-line -1)
      (skip-chars-forward " \t")
      (cond ((looking-at "\\*[ \t]+/")
	     ;; Delete region between `*' and `/' if there is only whitespaces.
	     (forward-char 1)
	     (delete-horizontal-space))
	    ((and (looking-at "/\\*\\|\\*[^/]")
		  (not (save-excursion
			 (search-forward "\\*/" (cocoa5-get-end-of-line) t))))
	     (setq setstar t))))
    ;; If last line was a star comment line then this one shall be too.
;    (if (null setstar)
;	(cocoa5-indent-line)
;      (insert "*  "))
))

(defun cocoa5-electric-semi-or-dot ()
  "Insert `;' or `.' character and, if in code, reindent the line."
  (interactive)
  (insert last-command-event)
  ;; Special handling if point is inside code
  (when (cocoa5-point-in-code-p)
    (save-excursion
      (if abbrev-mode (expand-abbrev))
      (beginning-of-line)
      (cocoa5-indent-line))
    (if cocoa5-auto-newline
	(cocoa5-electric-terminate-line))))

;; JAA commented out 2014-03-23
;; (defun cocoa5-electric-semi-or-dot ()
;;   "Insert `;' or `.' character and reindent the line."
;;   (interactive)
;;   (insert last-command-event)
;;   (if abbrev-mode (expand-abbrev))
;;   ;; Do nothing if within string.
;;   (if (cocoa5-within-string)
;;       ()
;;     (save-excursion
;;       (beginning-of-line)
;;       (cocoa5-indent-line))
;;     (if cocoa5-auto-newline
;; 	(cocoa5-electric-terminate-line))))

;; (defun cocoa5-electric-colon ()
;;   "Insert `:' and do all indentions except line indent on this line."
;;   (interactive)
;;   (insert last-command-event)
;;   ;; Do nothing if within string.
;;   (if (cocoa5-within-string)
;;       ()
;;     (save-excursion
;;       (beginning-of-line)
;;       (cocoa5-indent-line))
;;     (let ((cocoa5-tab-always-indent nil))
;;       (cocoa5-indent-command))))

(defun cocoa5-electric-equal ()
  "Insert `=', and do indention if within type declaration."
  (interactive)
  (insert last-command-event)
  (if (eq (car (cocoa5-calculate-indent)) 'declaration)
      (let ((cocoa5-tab-always-indent nil))
	(cocoa5-indent-command))))

;; (defun cocoa5-electric-hash ()
;;   "Insert `#', and indent to column 0 if this is a CPP directive."
;;   (interactive)
;;   (insert last-command-event)
;;   (if (save-excursion (beginning-of-line) (looking-at "^[ \t]*#"))
;;       (save-excursion (beginning-of-line)
;; 		      (delete-horizontal-space))))

(defun cocoa5-electric-tab-no-completion ()
  "Function called when TAB is pressed in CoCoA mode."
  (interactive)
  ;; Do nothing if within a string or in a CPP directive.
  (if (or (cocoa5-within-string)
	  (and (not (bolp))
	       (save-excursion (beginning-of-line) (eq (following-char) ?#))))
      (insert "\t")
    ;; If cocoa5-tab-always-indent, indent the beginning of the line.
    (if cocoa5-tab-always-indent
	(save-excursion
	  (beginning-of-line)
	  (cocoa5-indent-line))
      (if (save-excursion
	    (skip-chars-backward " \t")
	    (bolp))
	  (cocoa5-indent-line)
	(insert "\t")))
    (cocoa5-indent-command))
)


(defun cocoa5-electric-tab ()
  (interactive)
  (if cocoa5-tab-is-completion
      (progn (cocoa5-electric-tab-no-completion) (dabbrev-expand()))
    (cocoa5-electric-tab-no-completion)
    )
  )



(defun cocoa5-beg-of-defun ()
  "Move backward to the beginning of the current function or procedure."
  (interactive)
  (while (not (looking-at cocoa5-defun-re))
    (cocoa5-last-unended-begin))
)

(defun cocoa5-end-of-statement ()
  "Move forward to end of current statement."
  (interactive)
  (let ((parse-sexp-ignore-comments t)
	(nest 0) pos
	(regexp (concat "\\(" cocoa5-beg-block-re "\\)\\|\\("
			cocoa5-end-block-re "\\)")))
    (if (not (looking-at "[ \t\n]")) (forward-sexp -1))
    (or (looking-at cocoa5-beg-block-re)
	;; Skip to end of statement
	(setq pos (catch 'found
		    (while t
		      (forward-sexp 1)
		      (cond ((looking-at "[ \t]*;")
			     (skip-chars-forward "^;")
			     (forward-char 1)
			     (throw 'found (point)))
			    ((save-excursion
			       (forward-sexp -1)
			       (looking-at cocoa5-beg-block-re))
			     (goto-char (match-beginning 0))
			     (throw 'found nil))
			    ((eobp)
			     (throw 'found (point))))))))
    (if (not pos)
	;; Skip a whole block
	(catch 'found
	  (while t
	    (re-search-forward regexp nil 'move)
	    (setq nest (if (match-end 1)
			   (1+ nest)
			 (1- nest)))
	    (cond ((eobp)
		   (throw 'found (point)))
		  ((= 0 nest)
		   (throw 'found (cocoa5-end-of-statement))))))
      pos)))

;;;
;;; Other functions
;;;
(defun cocoa5-set-auto-comments ()
  "Insert `{ case }' or `{ NAME }' on this line if appropriate.
Insert `{ case }' if there is an `end' on the line which
ends a case block.  Insert `{ NAME }' if there is an `end'
on the line which ends a function or procedure named NAME."
  (save-excursion
    (forward-line -1)
    (skip-chars-forward " \t")
    (if (and (or (looking-at "\\<EndDefine;") (looking-at "\\<End;") )
	     (not (save-excursion
		    (end-of-line)
		    (search-backward "--" (cocoa5-get-beg-of-line) t))))
	(let ((type (car (cocoa5-calculate-indent))))
	  (if (eq type 'declaration)
	      ()
	      (let ((found-define nil))
		;; Check if this is the end of a function
		(save-excursion
		  (cocoa5-beg-of-defun)
		  (setq define-point (point))
		  (setq found-define (looking-at cocoa5-defun-re)))
		(if found-define
		    (progn
		      (end-of-line)
		      (delete-horizontal-space)
		      (insert " -- ")
		      (let (b e)
			(save-excursion
			  (setq b (progn (goto-char define-point)
					 (skip-chars-forward "^[:space:]") ;; skip over "define" keyword ;; skipping "a-zA-Z" did not work?!?
					 (skip-chars-forward "[:space:]")  ;; skip white space between "define" and fn name
					 (point))
				e (progn (skip-chars-forward "a-zA-Z0-9_") ;; skip over fn name
					 (point))))
			(insert-buffer-substring (current-buffer) b e))
		      (insert "")))))))))



;;;
;;; Indentation
;;;
(defconst cocoa5-indent-alist
  '((block . (+ ind cocoa5-indent-level))
    (case . (+ ind cocoa5-case-indent))
    (caseblock . ind) (cpp . 0)
    (declaration . (+ ind cocoa5-indent-level))
    (paramlist . (cocoa5-indent-paramlist t))
    (comment . (cocoa5-indent-comment t))
    (defun . ind) (contexp . ind)
    (unknown . 0) (string . 0)))

(defun cocoa5-indent-command ()
  "Indent for special part of code."
  (let* ((indent-str (cocoa5-calculate-indent))
	 (type (car indent-str))
	 (ind (car (cdr indent-str))))
    (cond ((and (eq type 'paramlist)
		(or (memq 'all cocoa5-auto-lineup)
		    (memq 'paramlist cocoa5-auto-lineup)))
	   (cocoa5-indent-paramlist)
	   (cocoa5-indent-paramlist))
	  ((and (eq type 'declaration)
		(or (memq 'all cocoa5-auto-lineup)
		    (memq 'declaration  cocoa5-auto-lineup)))
	   (cocoa5-indent-declaration))
	  ((and (eq type 'case) (not (looking-at "^[ \t]*$"))
		(or (memq 'all cocoa5-auto-lineup)
		    (memq 'case cocoa5-auto-lineup)))
	   (cocoa5-indent-case)))
    (if (looking-at "[ \t]+$")
	(skip-chars-forward " \t"))))

(defun cocoa5-indent-line ()
  "Indent current line as a CoCoA statement."
  (let* ((indent-str (cocoa5-calculate-indent))
	 (type (car indent-str))
	 (ind (car (cdr indent-str))))
    (if (looking-at "^[0-9a-zA-Z]+[ \t]*:[^:=]")
	(search-forward ":" nil t))
    (delete-horizontal-space)
    ;; Some things should not be indented
    (if (or (and (eq type 'declaration) (looking-at cocoa5-declaration-re))
	    (eq type 'cpp)
	    (looking-at cocoa5-defun-re))
	()
      ;; Other things should have no extra indent
      (if (looking-at cocoa5-noindent-re)
	  (indent-to ind)
	;; But most lines are treated this way:
	(indent-to (eval (cdr (assoc type cocoa5-indent-alist))))
	))))

(defun cocoa5-calculate-indent ()
  "Calculate the indent of the current CoCoA line.
Return a list of two elements: (INDENT-TYPE INDENT-LEVEL)."
  (save-excursion
    (let* ((parse-sexp-ignore-comments t)
	   (oldpos (point))
	   (state (save-excursion (parse-partial-sexp (point-min) (point))))
	   (nest 0) (par 0) (complete (looking-at "[ \t]*[Ee]nd\\>"))
	   (elsed (looking-at "[ \t]*[Ee]lse\\>"))
	   (type (catch 'nesting
		   ;; Check if inside a string, comment or parenthesis
		   (cond ((nth 3 state) (throw 'nesting 'string))
			 ((nth 4 state) (throw 'nesting 'comment))
			 ((> (car state) 0)
;;			  (goto-char (scan-lists (point) -1 (car state)))
			  (backward-up-list)
			  (setq par (1+ (current-column)))
			  (throw 'nesting 'contexp)
			  )
			 ((save-excursion (beginning-of-line)
					  (eq (following-char) ?#))
			  (throw 'nesting 'cpp)))
		   ;; Loop until correct indent is found
		   (while t
		     (if (bobp) (throw 'nesting 'unknown))
		     (backward-sexp 1)
		     (cond (;--Escape from case statements
			    (and (looking-at "[A-Za-z0-9]+[ \t]*:[^:=]")
				 (not complete)
				 (save-excursion (skip-chars-backward " \t")
						 (bolp))
				 (= (save-excursion
				      (end-of-line) (backward-sexp) (point))
				    (point))
				 (> (save-excursion (goto-char oldpos)
						    (beginning-of-line)
						    (point))
				    (point)))
			    (throw 'nesting 'caseblock))
			   (;--Nest block outwards
			    (looking-at cocoa5-beg-block-re)
			    (if (= nest 0)
				(cond ((looking-at "NULLA_DA_FISSARE\\>")
				       (throw 'nesting 'case))
				      ((looking-at "NULLA_DA_FISSARE\\>")
				       (throw 'nesting 'declaration))
				      (t (throw 'nesting 'block)))
			      (setq nest (1- nest))))
			   (;--Nest block inwards
			    (looking-at cocoa5-end-block-re)
			    (if (and (looking-at "[Ee]nd\\s ")
				     elsed (not complete))
				(throw 'nesting 'block))
			    (setq complete t
				  nest (1+ nest)))
			   (;--Defun (or parameter list)
			    (looking-at cocoa5-defun-re)
			    (if (= 0 par)
				(throw 'nesting 'defun)
			      (setq par 0)
			      (let ((n 0))
				(while (re-search-forward

	"\\(\\<NULLA_DA_FISSARE\\>\\)\\|\\<NULLA_DA_FISSARE\\>"
					oldpos t)
				  (if (match-end 1)
				      (setq n (1+ n)) (setq n (1- n))))
				(if (> n 0)
				    (throw 'nesting 'declaration)
				  (throw 'nesting 'paramlist)))))
			   (;--Declaration part
			    (looking-at cocoa5-declaration-re)
			    (if (save-excursion
				  (goto-char oldpos)
				  (forward-line -1)
				  (looking-at "^[ \t]*$"))
				(throw 'nesting 'unknown)
			      (throw 'nesting 'declaration)))
			   (;--If, else or while statement
			    (and (not complete)
				 (looking-at cocoa5-subblock-re))
			    (throw 'nesting 'block))
			   (;--Found complete statement
			    (save-excursion (forward-sexp 1)
					    (= (following-char) ?\;))
			    (setq complete t))
			   (;--No known statements
			    (bobp)
			    (throw 'nesting 'unknown))
			   )))))

      ;; Return type of block and indent level.
      (if (> par 0)                               ; Unclosed Parenthesis
	  (list 'contexp par)
	(list type (cocoa5-indent-level))))))

(defun cocoa5-indent-level ()
  "Return the indent-level the current statement has.
Do not count labels, case-statements or records."
  (save-excursion
    (beginning-of-line)
    (if (looking-at "[ \t]*[0-9a-zA-Z]+[ \t]*:[^:=]")
	(search-forward ":" nil t)
      (if (looking-at ".*=[ \t]*NULLA_DA_FISSARE\\>")
	  (search-forward "=" nil t)))
    (skip-chars-forward " \t")
    (current-column)))

(defun cocoa5-indent-comment (&optional arg)
  "Indent current line as comment.
If optional arg is non-nil, just return the
column number the line should be indented to."
;;     (let* ((stcol (save-excursion
;;		    (re-search-backward "--" nil t)
;; 		    (1+ (current-column)))))
;;       (if arg stcol
;; 	(delete-horizontal-space)
;; 	(indent-to stcol)))
;;  AMB 2014-08: it should NOT indent, but I don't know how to (not) do it!!!
  (indent-to (current-column)) ;; indent at first column
)


;from b to e nicely. The lineup string is str."

(defun cocoa5-get-lineup-indent (b e str)
  (save-excursion
    (let ((ind 0)
	  (reg (concat str "\\|\\(\\<NULLA_DA_FISSARE\\>\\)"))
	  nest)
      (goto-char b)
      ;; Get rightmost position
      (while (< (point) e)
	(setq nest 1)
	(if (re-search-forward reg (min e (cocoa5-get-end-of-line 2)) 'move)
	    (progn
	      ;; Skip record blocks
	      (if (match-beginning 1)
		  (cocoa5-declaration-end)
		(progn
		  (goto-char (match-beginning 0))
		  (skip-chars-backward " \t")
		  (if (> (current-column) ind)
		      (setq ind (current-column)))
		  (goto-char (match-end 0)))))))
      ;; In case no lineup was found
      (if (> ind 0)
	  (1+ ind)
	;; No lineup-string found
	(goto-char b)
	(end-of-line)
	(skip-chars-backward " \t")
	(1+ (current-column))))))

;;-------------------------------------------------------------->
;; cocoa5-close-block  {like tex-close-latex-block} C-c C-e
;;-------------------------------------------------------------->

;; (defun cocoa5-last-unended-begin ()
;;   "Leave point at the beginning of the last command that is unended."
;; 					;  (while (and (re-search-backward "\\(\\\\begin\\s *{\\)\\|\\(\\\\end\\s *{\\)")
;; 					;              (looking-at "\\\\end{"))
;;   (save-restriction
;;     (narrow-to-region (point-min) (point))
;;     (while (and (re-search-backward 
;; 		 (concat 
;; 		  "\\<\\(" cocoa5-startkeywords-list "\\|"
;; 		  cocoa5-endkeywords-list "\\)\\>"
;; 		  ;;		"[^\"]*\\(\"[^\"]*\"[^\"]*\\)*\\($\\|\\'\\)" ;; this just counts double-quote chars
;; 		  "[^\"]*\\(\"\\([^\\\"]\\|\\\\.\\)*\"[^\"]*\\)*\\($\\|\\'\\)" ;; skips strings (inside a string if there is a backslash it skips both the backslash and the next char).
;; 		  )
;; 		 )
;; 		(looking-at cocoa5-endkeywords-list))
;;       (cocoa5-last-unended-begin))))

(defun cocoa5-last-unended-begin ()
  "Leave point at the beginning of the last command that is unended."
  (interactive)
  (while (and
	  (cocoa5-last-start-or-end-keyword)
	  (looking-at cocoa5-end-block-re))
    (cocoa5-last-unended-begin)))

(defun cocoa5-close-block ()
  "Creates an End... to match the last unclosed command."
  (interactive "*")
  (let ((new-line-needed (bolp))
	text indentation)
    (save-excursion
      (condition-case nil
          (cocoa5-last-unended-begin)
        (error (error "Couldn't find unended command"))) ;; yes, (error (error ...)) is RIGHT!!
      (setq indentation (current-column))
;;      (re-search-forward "[ \n]*\\(\\s *[^)\n ]*\\)")
      (re-search-forward "\\([A-Za-z]+\\)")
      (setq text (buffer-substring (match-beginning 1) (match-end 1))))
    (indent-to indentation)
    (if (string< "a" text)  (insert "end" text ";")  (insert "End" text ";"))
    (cocoa5-electric-terminate-line)))
;    (if new-line-needed (insert ?\n))))


(defun cocoa5-last-start-or-end-keyword ()
  "Move point to start of last keyword (outside string/comment) before current point"
  (let ((keyword-regexp (concat 
			 "\\<\\(" cocoa5-startkeywords-list "\\|"
			 cocoa5-end-block-re "\\)\\>"
			 )))
    (while (and 
	    (re-search-backward keyword-regexp)
	    (not (cocoa5-point-in-code-p)))
      nil)
    t ;; final outcome, must be non-nil
    ))


;; (defun cocoa5-point-in-code-p ()
;;   "prototype"
;;   (interactive)
;;   (setq string-literal "\"\\([^\\\"]\\|\\\\.\\)*\"")
;;   (setq no-comment-or-string "\\([^\"/-]\\|/[^/]\\|-[^-]\\)*")
;;   (setq in-code (concat "^" no-comment-or-string "\\(" string-literal no-comment-or-string "\\)*" "\\'"))
;;   (save-restriction
;;     (narrow-to-region (point-min) (point))
;;     (save-excursion
;;       (beginning-of-line)
;;       ;;(re-search-forward in-code)
;;       (re-search-forward in-code (point-max) t)
;;       )
;;     )
;;   )

;; Messy: assumes start of curr line is not in multiline comment/string
;;        checks whether point is inside a string or a comment (or after '?')
(defun cocoa5-point-in-code-p ()
  "say whether point is not in a string or comment"
  ;;  (interactive)
  (progn
    (setq inline-comment "/\\*\\([^*]\\|\\*[^/]\\)*\\*/") ;; matches /*...stuff...*/
    (setq end-of-line-comment "\\(//\\|--\\)") ;; matches // or --
    (setq string-literal "\"\\([^\\\"]\\|\\\\.\\)*\"") ;; matches "...string..."
    (setq string-or-comment (concat "\\(" string-literal "\\|" inline-comment "\\)"))
    (setq no-string-or-comment "\\([^-\"/?]\\|/[^/*]\\|-[^-]\\)*") ;; matches until ?, ", /*, //, or --
    ;;  (setq online-help "\s-*\?") ;; \s- matches whitespace
    (setq in-code (concat "^" no-string-or-comment "\\(" string-or-comment no-string-or-comment "\\)*" "\\'"))
    (save-restriction
      (narrow-to-region (point-min) (point))
      (save-excursion
	(beginning-of-line)
	(setq last-point-in-code (re-search-forward in-code (point-max) t))
	(eq last-point-in-code (point-max))
))))


;;-------------------------------------------------------------;;
;;                        RUNNING CoCoA                        ;;
;;-------------------------------------------------------------;;

(defun cocoa5 ()
  "CoCoA interface.
The function splits the current window into two parts, one containing
a CoCoA shell and the other containing a file chosen by the user
in cocoa5 mode.

The file in cocoa5 mode can be modified like any text file; in addition
to that, there is are special key bindings (C-c ...) for sending a line,
a region or the whole buffer as input to CoCoA.  This saves the effort
of cutting and pasting lines of input.

Type  C-h f cocoa5-mode  to know about the other functions available in
cocoa5 mode.

Customization:
- change the value of the global variable cocoa5-dir to the 
  directory where your CoCoA input files are;
- change the value of cocoa5-shell-launch-dir to the directory where the 
  CoCoA executable is.
"
;  (setq default-directory cocoa5-shell-launch-dir)
  (interactive)
  (require 'comint)
  (require 'easymenu)  ;; MENU

;;;; uncomment the following lines if you have problems with trailing ^M
;;;; you can also add them to your .emacs file
;;(if (null system-uses-terminfo) nil
;;    (autoload 'shell-strip-ctrl-m "shell" nil t)
;;    (add-hook 'comint-output-filter-functions 'shell-strip-ctrl-m)
;;    (add-hook 'comint-input-filter-functions 'shell-ctrl-m)
;;)
  (cocoa5-make-shell)
  ;;  (delete-other-windows)
  ;;  (save-excursion
  (switch-to-buffer (process-buffer (get-process "cocoa5")))
  (goto-char (point-max))
  ;;  )
)

;; (defun cocoa5server ()
;; ;  (setq default-directory cocoa5-shell-launch-dir)
;;   (interactive)
;;   (require 'comint)
;;   (require 'easymenu)  ;; MENU

;; ;;;; uncomment the following lines if you have problems with trailing ^M
;; ;;;; you can also add them to your .emacs file
;; ;;(if (null system-uses-terminfo) nil
;; ;;    (autoload 'shell-strip-ctrl-m "shell" nil t)
;; ;;    (add-hook 'comint-output-filter-functions 'shell-strip-ctrl-m)
;; ;;    (add-hook 'comint-input-filter-functions 'shell-ctrl-m)
;; ;;)
;;   (cocoa5-make-server-shell)
;;   (switch-to-buffer-other-frame "*CoCoA5Server*")
;;   (previous-multiframe-window) ;possibly (other-frame -1)
;; )

;;----------------------------------------------------------------------
;; cocoa5-mode-map
;;----------------------------------------------------------------------
(message "cocoa5-mode-map")

(if (and (boundp 'cua-mode) cua-mode)
    (if (cua-mode)  ;; <------- ugly trick?  why do I need this line?????
	(progn
	  (message "cocoa5: cua-mode is on")
	  (setq cua-rectangle-mark-key (kbd "C-S-<return>"))
	  (cua-mode 1)
	  )
      )
  (message "cocoa5: cua-mode is off")
  )

(defvar cocoa5-mode-map ()  "Keymap used in CoCoA mode.")
(if cocoa5-mode-map
    ()
;; syntax
  (setq cocoa5-mode-map (make-sparse-keymap))
  (define-key cocoa5-mode-map ";"       'cocoa5-electric-semi-or-dot)
;;  (define-key cocoa5-mode-map "."       'cocoa5-electric-semi-or-dot)
  (define-key cocoa5-mode-map "="       'cocoa5-electric-equal)
  (define-key cocoa5-mode-map "\r"      'cocoa5-electric-terminate-line)
  (define-key cocoa5-mode-map "\t"      'cocoa5-electric-tab)
  (define-key cocoa5-mode-map "\177"    'backward-delete-char-untabify)
  (define-key cocoa5-mode-map "\M-\C-a" 'cocoa5-beg-of-defun)
;; to send commands to CoCoA
  (define-key cocoa5-mode-map "\C-c\C-l" 'cocoa5-send-line)
  (define-key cocoa5-mode-map (kbd "C-c RET") 'cocoa5-send-line)
  (define-key cocoa5-mode-map (kbd "C-<return>") 'cocoa5-send-line-or-region)
  (define-key cocoa5-mode-map (kbd "M-<return>") 'cocoa5-send-line-or-region)
  (define-key cocoa5-mode-map "\C-c\C-r" 'cocoa5-send-region)
  (define-key cocoa5-mode-map "\C-c\C-f" 'cocoa5-source-file)
  (define-key cocoa5-mode-map "\C-c\C-p" 'cocoa5-pop-to-buffer-source-find-error)
;; other useful things
  (define-key cocoa5-mode-map "\C-c\C-c" 'cocoa5-interrupt-or-comment-region)
  (define-key cocoa5-mode-map "\C-c\C-u" 'uncomment-region)
  (define-key cocoa5-mode-map "\C-c\C-e" 'cocoa5-close-block)
  (define-key cocoa5-mode-map "\C-c ]"   'cocoa5-close-block)
  (define-key cocoa5-mode-map "\C-c\C-o" 'cocoa5-toggle-selective-display-col)
  (define-key cocoa5-mode-map "\C-c\C-h" 'cocoa5-word-help)
;  (define-key cocoa5-mode-map "\C-cm"    'cocoa5-word-help)
;  (define-key cocoa5-mode-map "\C-c?"    'cocoa5-word-help)
;  (define-key cocoa5-mode-map "\C-c\C-n" 'cocoa5-open-file)
;  (define-key cocoa5-mode-map "\C-c\C-o" 'cocoa5-shell-panic-string)
  (define-key cocoa5-mode-map "\C-c\C-q" 'cocoa5-shell-quit)
  (define-key cocoa5-mode-map "\C-c\C-k" 'cocoa5-restart-shell)
;  (define-key cocoa5-mode-map "\C-c\C-w" 'cocoa5-open-wordlist-file)
;;  (define-key cocoa5-mode-map "\C-c " 'cocoa5-rotate-split)
;;  (define-key cocoa5-mode-map (kbd "C-c C-SPC") 'cocoa5-rotate-split)
  (define-key cocoa5-mode-map "\C-c\C-s" 'cocoa5-side-by-side)
  (define-key cocoa5-mode-map "\C-c\C-t" 'cocoa5-top-and-bottom)
)

;;-------------------------------------------------------------->
;; cocoa5-menu
;;-------------------------------------------------------------->
(message "cocoa5-menu")

(defvar cocoa5-mode-menu nil)
(defvar cocoa5-menu nil  "Menu for `cocoa5-mode'.")
(setq cocoa5-menu
      '("CoCoA-5"
	["Send Line to CoCoA-5"        cocoa5-send-line   t]
	["Source Region to CoCoA-5"    cocoa5-send-region t]
	["Source File into CoCoA-5"    cocoa5-source-file t]
	["Help for Word"               cocoa5-word-help]
	["Go to Error"                 cocoa5-pop-to-buffer-source-find-error t]
	["CoCoA\&source: side-by-side"    cocoa5-side-by-side
	 :style toggle
	 :selected cocoa5-split-is-side-by-side
	 :active 't]
	["CoCoA\&source: top-and-bottom"  cocoa5-top-and-bottom
	 :style toggle
	 :selected (not cocoa5-split-is-side-by-side)
	 :active 't]
	["Enable SystemCommand on start-up"    cocoa5-enable-system-command
	 :style toggle
	 :selected cocoa5-enable-system-cmd-flag
	 :active 't]
	"---"
	["Indent Region"                  indent-region]
	["Comment Region"                 comment-region]
	["Uncomment Region"               uncomment-region]
	["Outline from this column"   (cocoa5-toggle-selective-display-col)
	 :style toggle
	 :selected selective-display
	 ]
;; 	["Outline functions"  cocoa5-toggle-selective-display
;; 	 :style toggle
;; 	 :selected selective-display
;; 	 ]
	"---"
	["Print EndIf/EndFor/.."          cocoa5-close-block]
	["Word Completion"                dabbrev-expand]
	"---"
	["Echo Source Region to CoCoA-5" (setq cocoa5-echo-source-region
					    (not cocoa5-echo-source-region))
	 :style toggle
	 :selected cocoa5-echo-source-region
	 :active 't]
	["Use TAB for Word Completion"      (setq cocoa5-tab-is-completion
					    (not cocoa5-tab-is-completion))
	 :style toggle
	 :selected cocoa5-tab-is-completion
	 :active 't]
	["Capitalize keywords"      abbrev-mode
	 :style toggle
	 :selected abbrev-mode
	 :active 't]
	["Preferences"                     cocoa5-customize]
	"---"
	["Show cocoa5 version"             cocoa5-show-version t]
	["Show cocoa5-mode version"        cocoa5-show-mode-version t]
	"---"
	["Interrupt CoCoA5"                cocoa5-interrupt-shell t]
	["(Re)start CoCoA5"                cocoa5-restart-shell t]
	["Quit CoCoA-5 and Emacs"          cocoa5-quit-emacs t]
	))

(message "cocoa5-comint-menu")

(defvar cocoa5-comint-mode-menu nil)
(defvar cocoa5-comint-menu nil  "Menu for `cocoa5-comint-mode'.")
(setq cocoa5-comint-menu
      '("CoCoA-5"
	["Help for Word"             cocoa5-word-help]
	["Go to Error"           cocoa5-pop-to-buffer-source-find-error t]
	"---"
	["Print EndIf/EndFor/.."          cocoa5-close-block]
	["Word Completion"                dabbrev-expand]
	"---"
	["Preferences"                     cocoa5-customize]
	["Show cocoa5-mode version"        cocoa5-show-mode-version t]
	["Show cocoa5 version"             cocoa5-show-version t]
	"---"
	["Interrupt CoCoA5"                cocoa5-interrupt-shell t]
	["(Re)start CoCoA5"                 cocoa5-comint-restart-shell t]
	["Enable SystemCommand on start-up"    cocoa5-enable-system-command
	 :style toggle
	 :selected cocoa5-enable-system-cmd-flag
	 :active 't]
	["Quit CoCoA-5 and Emacs"           cocoa5-quit-emacs t]
	))

;; (defvar cocoa5-statements-menu nil  "Menu for cocoa5 statements.")
;; (setq cocoa5-statements-menu
;;       '("Statements"
;; 	["If"   cocoa5-if-statement]
;; 	["For"   cocoa5-for-statement]
;; 	["Foreach"   cocoa5-foreach-statement]
;; 	["While"   cocoa5-while-statement]
;; ))


;; (defun cocoa5-while-statement ()
;; ""
;; (insert
;; "While . Do
;;   .
;; EndWhile"
;; )
;; )

(defun cocoa5-mode-setup-menubar ()
  "Initial setup of cocoa5 and insertions menus."
  (progn
    (easy-menu-define			; set up cocoa5 menu
      cocoa5-mode-menu cocoa5-mode-map "Menu used in cocoa5-mode"
      cocoa5-menu)
    (easy-menu-add cocoa5-mode-menu cocoa5-mode-map) )
)

(defun cocoa5-comint-setup-menubar ()
  "Initial setup of cocoa5 and insertions menus."
  (progn
    (easy-menu-define			; set up menu
      cocoa5-comint-mode-menu cocoa5-comint-mode-map "Menu used in cocoa5-comint-mode"
      cocoa5-comint-menu)
    (easy-menu-add cocoa5-comint-mode-menu cocoa5-comint-mode-map) )
)



(defconst cocoa5-xemacs-p (string-match "XEmacs" (emacs-version)))

;;-------------------------------------------------------------->
;; cocoa-outline
;;-------------------------------------------------------------->
;; (defun cocoa5-outline ()
;;   "The function opens a file outline"
;;   (interactive)
;;   (setq cocoa-source-file buffer-file-name)
;;   (grep  (concat "grep -nH -e \"\\<Define\\>\" " cocoa5-source-file))
;; )

;;-------------------------------------------------------------->
;; cocoa-customize
;;-------------------------------------------------------------->

(defun cocoa5-customize ()
  "The function opens customize-group cocoa5"
  (interactive)
  (customize-group 'cocoa5)
)

;;-------------------------------------------------------------->
;; cocoa5-open-wordlist-file
;;-------------------------------------------------------------->

(defun cocoa5-open-wordlist-file ()
  "The function opens a wordlist file chosen by the user in cocoa5-mode"
  (interactive)
  (setq cocoa5-wordlist-file 
	(read-file-name "CoCoA-5 word list file: " cocoa5-wordlist-default-file))
  (if cocoa5-wordlist-file  
      (setq cocoa5-wordlist-file cocoa5-wordlist-default-file))
  (find-file-noselect cocoa5-wordlist-file)
)

;;-------------------------------------------------------------->
;; cocoa5-open-file
;;-------------------------------------------------------------->

(setq cocoa5-shell-buffer-list '())

;; (defun cocoa5-open-file ()
;;   "The function opens a new file chosen by the user in cocoa5-mode"
;;   (interactive)
;;   (setq cocoa5-file 
;; 	(read-file-name "Enter name of a CoCoA-5 file: " cocoa5-dir))
;;   (setq cocoa5-file-buffer-name (find-file-noselect cocoa5-file))
;;   (switch-to-buffer cocoa5-file-buffer-name)
;;   (cocoa5-mode)
;;   (setq cocoa5-shell-buffer-list 
;; 	(cons cocoa5-file-buffer-name cocoa5-shell-buffer-list))
;; )

;;-------------------------------------------------------------->
;; cocoa5-make-shell
;;-------------------------------------------------------------->
(message "cocoa5-make-shell")

(defun cocoa5-make-shell ()
  (save-excursion
    ;;    (setq default-directory cocoa5-shell-launch-dir)

    (if (not (file-exists-p cocoa5-executable))
	(error (concat "Cannot find CoCoA-5 executable: " cocoa5-executable)))
    (if (not (file-executable-p cocoa5-executable))
	(error (concat "CoCoA-5 exists but is not executable: " cocoa5-executable)))

    (setq cocoa5-args  
	  (list
	   cocoa5-warning-level
	   "--no-readline"
	   "--prompt" cocoa5-prompt) )
    (if cocoa5-enable-system-cmd-flag
	(setq cocoa5-args (cons "--enable-system-cmd" cocoa5-args)))
    (if cocoa5-full-cocoalib-error
	(setq cocoa5-args (cons "--fullCoCoALibError" cocoa5-args)))
    (if (eq system-type 'windows-nt) 
	(setq cocoa5-args (cons "--packageDir"
				(cons 
				 (concat cocoa5-emacs-dir "../packages") 
				 cocoa5-args) )))

    (set-buffer  (apply 'make-comint
			"cocoa5" cocoa5-executable nil
			cocoa5-args))
   (if (not (eq (buffer-local-value 'major-mode (get-buffer "*cocoa5*")) `cocoa5-comint-mode))
       (cocoa5-comint-mode))
    (sleep-for 0.35)	 ;; allow CoCoA-5 some time to get started
    ;; Load wordlist
    (if cocoa5-auto-load-wordlist
	(find-file-noselect cocoa5-wordlist-default-file 'NOWARN))
    )
  )

;;-------------------------------------------------------------->
;; cocoa5-mode
;;-------------------------------------------------------------->

;;;###autoload
(defun cocoa5-mode ()
  "Major mode for editing a cocoa5 file that sends input to CoCoA.
Editing mode with \\[cocoa5-send-line] sending the current line as input.

Special commands:
\\{cocoa5-mode-map}
"
;;  (setq default-directory cocoa5-shell-launch-dir)
  (interactive)
  (kill-all-local-variables)
  (use-local-map cocoa5-mode-map)
  (setq major-mode 'cocoa5-mode
	mode-name "CoCoA5")

  (run-hooks 'cocoa5-mode-hook)
)

;;-------------------------------------------------------------->
;; cocoa5-send-string
;;-------------------------------------------------------------->

(defun cocoa5-send-string (cocoa5-input-line)
  "Send string cocoa5-input-line as input to CoCoA in the other window.
Note: if the *cocoa5* buffer has no process, CoCoA will be restarted
automatically when this function is called.
"
  (interactive)
  (cocoa5-shell-check-process)
;;  (end-of-line)
;;  (setq cocoa5-line-end (point))
;;  (setq cocoa5-last-line-without-newline (eobp))
;;  (beginning-of-line)
;;  (setq cocoa5-line-beg (point))
;;  (setq cocoa5-line-length (- cocoa5-line-end cocoa5-line-beg))
;;  (setq cocoa5-line-blankp (= (skip-chars-forward " \t") cocoa5-line-length))
;;  (if (not cocoa5-line-blankp)
;;      (if (> cocoa5-line-length 999) ;; use cocoa5-send-region for lines over 999 chars
;;					;then
;;	  (let ((save-echo cocoa5-echo-source-region))
;;	    (setq cocoa5-echo-source-region nil)
;;	    (cocoa5-send-region cocoa5-line-beg cocoa5-line-end)
;;	    (setq cocoa5-echo-source-region save-echo))
					;else
	(save-excursion
;;	  (setq cocoa5-input-line (buffer-substring cocoa5-line-beg cocoa5-line-end))
	  (setq cocoa5-source-buffer (current-buffer)) 
	  (cocoa5-shell-check-window)
	  (cocoa5-pop-to-buffer-comint)
	  (goto-char (point-max))
	  (move-beginning-of-line nil) (if (not (eq (point) (point-max))) (kill-line))
	  (insert cocoa5-input-line)
	  (comint-send-input)
	  (pop-to-buffer cocoa5-source-buffer)
	  ))

;;-------------------------------------------------------------->
;; cocoa5-send-line
;;-------------------------------------------------------------->

(defun cocoa5-send-line ()
  "Send current line as input to CoCoA in the other window.
Note that the cursor can be at any position in the line;
after the input has been sent, the cursor goes to the next non-empty
line in the cocoa file.

Note: if the *cocoa5* buffer has no process, CoCoA will be restarted
automatically when this function is called.
"
  (interactive)
  (cocoa5-shell-check-process)
  (end-of-line)
  (setq cocoa5-line-end (point))
  (setq cocoa5-last-line-without-newline (eobp))
  (beginning-of-line)
  (setq cocoa5-line-beg (point))
  (setq cocoa5-line-length (- cocoa5-line-end cocoa5-line-beg))
  (setq cocoa5-line-blankp (= (skip-chars-forward " \t") cocoa5-line-length))
  (if (not cocoa5-line-blankp)
      (let ((use-source-region (or (> cocoa5-line-length 999)   ;; check length first!
				   (string-match "[^[:alnum:][:space:][:punct:]]"  ;; REGEXP matches "funny chars" 
;;				   (string-match "[---ÿ]" ;; REGEXP matches unprintable ASCII (except TAB char)
						 (buffer-substring-no-properties cocoa5-line-beg cocoa5-line-end)))))

	(if use-source-region ;; if longer than 999 or contains funny chars
	    (let ((save-echo cocoa5-echo-source-region))
	      (setq cocoa5-echo-source-region nil)
	      (cocoa5-send-region cocoa5-line-beg cocoa5-line-end)
	      (setq cocoa5-echo-source-region save-echo))
					;else
	(save-excursion
	  (setq cocoa5-input-line (buffer-substring cocoa5-line-beg cocoa5-line-end))
	  (setq cocoa5-source-buffer (current-buffer)) 
	  (cocoa5-shell-check-window)
	  (cocoa5-pop-to-buffer-comint)
	  (goto-char (point-max))
	  (move-beginning-of-line nil) (if (not (eq (point) (point-max))) (kill-line))
	  (insert cocoa5-input-line)
	  (comint-send-input)
	  (pop-to-buffer cocoa5-source-buffer)
	  ))))
  (end-of-line)
  ;; if last line was non-empty and without newline, insert a newline
  (if (and cocoa5-last-line-without-newline
	   (not cocoa5-line-blankp))
      (insert "\n"))  ;; Emacs will report error if the buffer is read-only
  (skip-chars-forward " \t\n")
)

;;-------------------------------------------------------------->
;; cocoa5-show-version
;;-------------------------------------------------------------->

(defun cocoa5-exec-this (arg)
  "Send arg to CoCoA5
"
  (setq cocoa5-source-buffer (current-buffer)) ;  for jumping back to it. Burki
  (save-excursion
    (cocoa5-shell-check-process)
    (cocoa5-shell-check-window)
;    (switch-to-buffer-other-window "*cocoa5*")
;    (pop-to-buffer "*cocoa5*")
    (cocoa5-pop-to-buffer-comint)
    (goto-char (point-max))
    (move-beginning-of-line nil)  (if (not (eq (point) (point-max))) (kill-line))
    (insert arg)
    (comint-send-input)
    (pop-to-buffer cocoa5-source-buffer)
    )
  (skip-chars-forward "-*\n-*")
)

(defun cocoa5-show-version ()
  (interactive)
  (cocoa5-exec-this "indent(VersionInfo());")
)

(defun cocoa5-show-mode-version ()
  (interactive)
  (cocoa5-exec-this (concat "-- " cocoa5-mode-version))
)

;; (defun cocoa5-shell-panic-string ()
;;   (interactive)
;;   (cocoa5-exec-this "\"'*/ # \"CoCoA \"+/*NOT*/\"ready\"; -- repeat to get out of pending statements")
;; )

;;-------------------------------------------------------------->
;; cocoa5-word-help
;;-------------------------------------------------------------->

(defun cocoa5-word-help ()
  "Send current word input to CoCoA 'Man' in the *cocoa5* window.
Note that the cursor can be at any position in the word;
after the input has been sent, the cursor goes to the next non-empty
line in the cocoa5 file.
Note: if the *cocoa5* buffer has no process, CoCoA will be restarted
automatically when this function is called.
"
  (interactive)
  (cocoa5-shell-check-process)
  (setq word-regexp "[a-zA-Z][a-zA-Z0-9_]*")
  (setq word-boundary "\\(^\\|[^a-zA-Z_]\\)")
  ;; Now mark the word and save to string.
  (forward-char)
  (re-search-backward (concat word-boundary word-regexp))
  (or (re-search-forward word-regexp (point-max) t)
      (error "No word found to check!"))
  (setq start (match-beginning 0)
	end (point)
	word (buffer-substring start end))
  (cocoa5-exec-this "--v-v-v--( Help for word )--v-v-v--")
  (cocoa5-exec-this (concat "? " word))
)
;;-------------------------------------------------------------->
;; cocoa5-shell-quit
;;-------------------------------------------------------------->
(message "cocoa5-shell-quit")

(defun cocoa5-shell-quit ()
  "Quit CoCoA after saving all notebook files."
  (interactive)
  (message "saving all notebook files . . .")
  (save-excursion
    (while (not (null cocoa5-shell-buffer-list))
      (setq cocoa5-shell-buffer-tosave (car cocoa5-shell-buffer-list))
      (if (buffer-name cocoa5-shell-buffer-tosave) 
	  (progn 
	    (switch-to-buffer cocoa5-shell-buffer-tosave)
	    (save-buffer)
	    (kill-buffer cocoa5-shell-buffer-tosave)))
      (setq cocoa5-shell-buffer-list (cdr cocoa5-shell-buffer-list))
      ))
  (cocoa5-shell-quit-cocoa5-shell)
  (message "  . . . bye!  ")
)

;;-------------------------------------------------------------->
;; cocoa5-kill-shell-process
;;-------------------------------------------------------------->

(defun cocoa5-kill-shell-process ()
  (if (not (eq (get-buffer-process "*cocoa5*") nil)) 
      (delete-process "*cocoa5*"))
;;      (process-kill-without-query (get-process "cocoa5")))
)

;;-------------------------------------------------------------->
;; cocoa5-interrupt-shell-process
;;-------------------------------------------------------------->

(defun cocoa5-interrupt-shell-process ()
  (if (not (eq (get-buffer-process "*cocoa5*") nil)) 
      (comint-interrupt-subjob))
)

;;-------------------------------------------------------------->
;; cocoa5-quit-emacs
;;-------------------------------------------------------------->

(defun cocoa5-quit-emacs ()
"Quit cocoa5 and then emacs."

(interactive)
;;(process-kill-without-query (get-process "cocoa5"))
(cocoa5-kill-shell-process)
;(cocoa5-kill-server-process)
(save-buffers-kill-emacs)
)

;;-------------------------------------------------------------->
;; cocoa5-restart-shell
;;-------------------------------------------------------------->
(message "cocoa5-restart-shell")

(defun cocoa5-comint-restart-shell ()
"(re)start CoCoA-5 process in the buffer named *cocoa5*."

(interactive)
(cocoa5-kill-shell-process)
(cocoa5-make-shell)
)

(defun cocoa5-restart-shell ()
"(Re)start CoCoA process in the buffer named *cocoa5*."

(interactive)
(cocoa5-kill-shell-process)
(cocoa5-make-shell)
(switch-to-buffer-other-window "*cocoa5*")
(goto-char (point-max))
(other-window 1)
)

;;-------------------------------------------------------------->
;; cocoa5-interrupt-shell
(message "cocoa5-interrupt-shell")
;;-------------------------------------------------------------->

(defun cocoa5-comint-interrupt-shell ()
"Interrupt CoCoA-5 process in the buffer named *cocoa5*."

(interactive)
(cocoa5-interrupt-shell-process)
(cocoa5-make-shell)
)

(defun cocoa5-interrupt-shell ()
"Interrupt CoCoA process in the buffer named *cocoa5*."

(interactive)
(switch-to-buffer-other-window "*cocoa5*")
(cocoa5-interrupt-shell-process)
(other-window 1)
)

;;-------------------------------------------------------------->
;; cocoa5-region-active-p
;;-------------------------------------------------------------->
(message "cocoa5-region-active-p")
;; for emacs-22 (use-region-p first defined in emacs 23)
(defun cocoa5-region-active-p ()
  "Return t if Transient Mark mode is enabled and the mark is active."
  (and transient-mark-mode mark-active))

;;-------------------------------------------------------------->
;; cocoa5-comint-interrupt
;;-------------------------------------------------------------->
(message "cocoa5-comint-interrupt")

(defun cocoa5-comint-interrupt ()
  "Send interrupt to CoCoA-5 process in the buffer named *cocoa5*."

  (interactive)
  (cocoa5-make-shell)
  (switch-to-buffer-other-window "*cocoa5*")
  (goto-char (point-max))
  (other-window 1)
  )

(defun cocoa5-interrupt-or-comment-region ()
  "If a region is selected, comment-region;
otherwise interrupts CoCoA-5 process.
"
  (interactive)
;;  (if (use-region-p) ;; DOES NOT WORK WITH OLD EMACS versions
(if (if (fboundp 'use-region-p) (use-region-p) (cocoa5-region-active-p))
;then
      (comment-region (region-beginning) (region-end))
;else
    (cocoa5-interrupt-shell))
)


;;-------------------------------------------------------------->
;; cocoa5-send-region
;;-------------------------------------------------------------->

;; (defun cocoa5-send-region ()
;;   "Send current region as input to CoCoA5.

;; Note: if the *cocoa5* buffer has no process, CoCoA will be restarted
;; automatically when this function is called.
;; "
;;   (interactive)
;;   (cocoa5-shell-check-process)
;;   (setq kill-ring '())
;;   (copy-region-as-kill (point) (mark))
;;   (cocoa5-shell-check-window)
;;   (save-excursion
;;     (switch-to-buffer-other-window "*cocoa5*")
;;     (if (not (eq (process-status "cocoa5") 'run)) (cocoa5-make-shell))
;;     (goto-char (point-max))
;;     (yank)
;;     (comint-send-input)
;; ;    (switch-to-buffer-other-window cocoa5-file-buffer-name)
;;     (other-window 1)
;;   )
;; )

;;-------------------------------------------------------------->
;; cocoa5-send-region
;;-------------------------------------------------------------->

(defun cocoa5-char-index-in-line (pos)
  (save-excursion
    (progn
      (goto-char pos)
      (setq end-pt (point))
      (forward-line 0) ;; move to beginning of line
      (setq start-pt (point))
      (- end-pt start-pt))))


;; for  send-region
(defun cocoa5-comment-line-replace (region-in)
  "replace ^ with -- in the string region-in"
  (with-temp-buffer
    (insert region-in)
    (goto-char (point-min))
    (while (search-forward-regexp "^" nil t) (replace-match "-- " nil t))
;;    (buffer-substring (point-min) (- (point-max) 1)))  ;; skip the last newline
    (buffer-substring (point-min) (point-max)))
)

;; Returns a CoCoA-5 string literal containing the name of the file associated to the curr buffer
;; Note that the return string already contains an initial and final double quote!
(defun cocoa5-string-file-name ()
  (setq escaped-file-name buffer-file-name)
  (setq escaped-file-name (replace-regexp-in-string "\\\\" "\\\\\\\\" escaped-file-name))   ; escape \
  (setq escaped-file-name (replace-regexp-in-string "\"" "\\\\\"" escaped-file-name))       ; escape "
  (setq escaped-file-name (replace-regexp-in-string "\n" "\\\\n" escaped-file-name)); escape newline
  (concat "\"" escaped-file-name "\"")
)

(defun cocoa5-send-region (region-beg region-end)
  "Send current region as input to CoCoA5.

Note: if the *cocoa5* buffer has no process, CoCoA will be restarted
automatically when this function is called.
"
  (interactive "r")
  (cocoa5-shell-check-process)
  (setq kill-ring '())
  (save-buffer)
  (setq cocoa5-shell-source-command
    (concat "SourceRegion "
	    (number-to-string (line-number-at-pos region-beg))
	    ","
	    (number-to-string (1+ (cocoa5-char-index-in-line region-beg)))
	    " to "
	    (number-to-string (line-number-at-pos region-end))
	    ","
	    (number-to-string (1+ (cocoa5-char-index-in-line region-end)))
	    " in " (cocoa5-string-file-name)
	    ";"))
  (if cocoa5-echo-source-region
      (setq cocoa5-sourceregion-region 
	    (cocoa5-comment-line-replace
	     (buffer-substring region-beg region-end)
	     )))
  (save-excursion
    (setq cocoa5-source-buffer (current-buffer)) ;  for jumping back to it. Burki
    (cocoa5-shell-check-window)
    (cocoa5-shell-check-process)
    (cocoa5-pop-to-buffer-comint)
    (goto-char (point-max))
    (move-beginning-of-line nil)  (if (not (eq (point) (point-max))) (kill-line))
    (insert cocoa5-shell-source-command)
    (if cocoa5-echo-source-region
	(insert "\n" cocoa5-sourceregion-region) ;; output region
      )
    (comint-send-input)
    (pop-to-buffer cocoa5-source-buffer)
  )
  (goto-char region-end)
  (deactivate-mark)
)

;;-------------------------------------------------------------->
;; cocoa5-send-line-or-region
;;-------------------------------------------------------------->

(defun cocoa5-send-line-or-region ()
  "If a region is selected, send it as input to CoCoA;
otherwise send the current line to CoCoA."
  (interactive)
;;  (if (use-region-p) ;; DOES NOT WORK WITH OLD EMACS versions
(if (if (fboundp 'use-region-p) (use-region-p) (cocoa5-region-active-p))
;then
      (cocoa5-send-region (region-beginning) (region-end))
;else
    (cocoa5-send-line))
)

;;-------------------------------------------------------------->
;; cocoa5-source-file
;;-------------------------------------------------------------->

(defun cocoa5-source-file ()
  "Saves current file (no prompt!) and send current file's as input to CoCoA5.

Note: if the *cocoa5* buffer has no process, CoCoA will be restarted
automatically when this function is called.
"
  (interactive)
  (cocoa5-shell-check-process)
  (save-buffer)
  (setq cocoa5-shell-source-command (concat "source " (cocoa5-string-file-name) ";"))
  ;;    (if (cocoa5-shell-running)
  ;;        (cocoa5-shell-kill-job)
  ;;      (cocoa5-make-shell)
  ;;     )
  (cocoa5-shell-check-window)
  (save-excursion
    (setq cocoa5-source-buffer (current-buffer)) ;  for jumping back to it. Burki
    (cocoa5-shell-check-window)
    (cocoa5-shell-check-process)
    (cocoa5-pop-to-buffer-comint)
;;???    (if (not (eq (process-status "cocoa5") 'run)) (cocoa5-make-shell))
    (goto-char (point-max))
    (move-beginning-of-line nil)  (if (not (eq (point) (point-max))) (kill-line))
    (insert cocoa5-shell-source-command)
    (comint-send-input)
    (pop-to-buffer cocoa5-source-buffer)
  )
)

;;-------------------------------------------------------------->
;; cocoa5-shell-check-process
;;-------------------------------------------------------------->
(message "cocoa5-shell-check-process")

(defun cocoa5-shell-check-process ()
  (if (eq (get-buffer-process "*cocoa5*") nil) (cocoa5-make-shell) )
  ;;    (if (not (eq (process-status "cocoa5") 'run)) (cocoa5-make-shell))
)

(defun cocoa5-split-horizontally ()
  "Place cocoa5 source and shell buffers left and right"
  (delete-other-windows)
  (split-window-horizontally)
  (switch-to-buffer-other-window "*cocoa5*")
  (other-window 1)
  )

(defun cocoa5-split-vertically ()
  "Place cocoa5 shell above and source below"
  (delete-other-windows)
  (split-window-vertically)
  (switch-to-buffer "*cocoa5*")
  (other-window -1)
  )

(defun cocoa5-side-by-side ()
  "Place cocoa5 source and shell buffers side by side, and set flag"
  (interactive)
  (setq cocoa5-split-is-side-by-side t)
  (cocoa5-shell-check-process)
  (cocoa5-split-horizontally)
  )

(defun cocoa5-top-and-bottom ()
  "Place cocoa5 shell above and source below, and set flag"
  (interactive)
  (setq cocoa5-split-is-side-by-side nil)
  (cocoa5-shell-check-process)
  (cocoa5-split-vertically)
  )

(defun cocoa5-enable-system-command ()
  "When starting cocoa5 enable SystemCommand"
  (interactive)
  (setq cocoa5-enable-system-cmd-flag (not cocoa5-enable-system-cmd-flag))
  )

(defun cocoa5-shell-check-window ()
  "If a *cocoa5* buffer is open use it, otherwise call cocoa5-split-***"
  (if (get-buffer-window "*cocoa5*" 0)
      ()
    (if cocoa5-split-is-side-by-side (cocoa5-split-horizontally)
      (cocoa5-split-vertically))
    )
  )

;; ;; Found this fn on an emacs forum
;; ;; Should add some extra checks (e.g. just 2 subwindows)
;; (defun cocoa5-rotate-split ()
;;   "Rotate window split from vertical to horizontal or vice versa"
;;   (interactive)
;;   (let ((root (car (window-tree))))
;;     (if (listp root)
;;         (let* ((w1 (nth 2 root))
;;                (w2 (nth 3 root))
;;                (b1 (window-buffer w1))
;;                (b2 (window-buffer w2)))
;;           (cond ((car root)             ; currently vertically split
;;                  (delete-window w2)
;;                  (set-window-buffer (split-window-horizontally) b2))
;;                 (t                      ; currently horizontally split
;;                  (delete-window w1)
;;                  (set-window-buffer (split-window-vertically) b1))))
;;       (message "Root window not split")))) 



;; ;; installation: Append this file to your cocoa5.el file.

;; ;;-------------------------------------------------------------->
;; ;; cocoa5-comint-mode
;; ;; =
;; ;; comint-mode + things for CoCoA5.
;; ;;
;; ;; things = F2 key, ...
;; ;;
;; ;; Burkhard Zimmermann, Friday May 13, 2005.
;; ;;-------------------------------------------------------------->


;; ;; tested with: emacs:              GNU emacs 21.3.1, 
;; ;;              cocoa5-mode-version: "21 December 2004",
;; ;;              Windows XP.

;; (message "Adding F2 to CoCoA mode.")
;; (set 'cocoa5-mode-version (concat cocoa5-mode-version " - [F2]"))

;; ;; F2 = send file to CoCoA5, similarly to C-c C-a.
;; (add-hook 'cocoa5-mode-hook 
;;    (lambda () (define-key cocoa5-mode-map [f2] 'cocoa5-send-file-xxx)
;;               (define-key cocoa5-mode-map [f1] 'cocoa5-word-help))) ; Burki.
;; ; (I put the to the define-key in the hook so that I don't need to modify cocoa5.el at all.)

;; ;; change to the function: cocoa5-make-server-shell
;; (defun cocoa5-make-server-shell ()
;;   (save-excursion
;;     (set-buffer  (make-comint "CoCoA5Server" 
;; 			      (concat cocoa5server-executable)))
;;     (cocoa5-comint-mode) ;; <-- added line.
;;   ))

;; ;;-------------------------------------------------------------->
;; ;; cocoa5-send-file-xxx
;; ;;
;; ;; Similar to cocoa5-send-file. differences:
;; ;;
;; ;; - remembers cocoa5-source-buffer. that's used for jumping back.
;; ;; - gives focus to the comint buffer.
;; ;; - keeps the windows layout. (uses pop-to-buffer. doesn't call cocoa5-shell-check-window)
;; ;;-------------------------------------------------------------->

;; (defun cocoa5-send-file-xxx ()
;;   "Sends the content of the buffer to CoCoA5. Similar to cocoa5-send-file."
;;   (interactive)
;;   (cocoa5-shell-check-process)
;;   (save-buffer) ;; <--- ;; to do: don't save autoamtically. save in a temporal file only. 
;;   (setq cocoa5-shell-input-file buffer-file-name)
;;   (setq cocoa5-shell-input-file-name (concat (concat "<< '" cocoa5-shell-input-file) "';"))
;;   ;; I don't like to have my window layout changed to a horizontal split by:
;;   ;; (cocoa5-shell-check-window)
;;   (save-excursion
;;     (set 'cocoa5-source-buffer (current-buffer)) ; <---  for jumping back to it. ; Burki
;;     (pop-to-buffer "*cocoa5*")
;;     (if (not (eq (process-status "cocoa5") 'run)) (cocoa5-make-shell))
;;     (goto-char (point-max))
;;     (insert cocoa5-shell-input-file-name)
;;     (comint-send-input)
;;   ))


;; ;;-------------------------------------------------------------->
;; ;; cocoa5-comint-mode
;; ;; =
;; ;; comint-mode + things for CoCoA5.
;; ;;-------------------------------------------------------------->

(message "cocoa5-comint-mode")
(defun cocoa5-comint-mode ()
  "Major mode for interacting with CoCoA5. Similar to `comint-mode'.
\\{cocoa5-comint-mode-map}"
  (interactive)
  ;; stolen from maplev.el:
  (comint-mode)
  (setq comint-get-old-input (function cocoa5-get-old-input))
  (setq ;;comint-prompt-regexp (concat "^\\(" cocoa5-prompt "\\)+ *")
        ;; comint-use-prompt-regexp t
        comint-eol-on-send t ;; ?
        major-mode 'cocoa5-comint-mode
        mode-name  "cocoa5-comint")
  (use-local-map cocoa5-comint-mode-map)
  (cocoa5-comint-setup-menubar)
  (run-hooks 'cocoa5-comint-mode-hook)
  )

;; shouldn't this go into the function cocoa5-comint-mode?
(defvar cocoa5-comint-mode-map () "Keymap used in cocoa5-comint mode.")
(require 'comint)
(unless cocoa5-comint-mode-map
  (let ((map (copy-keymap comint-mode-map) ))
    (define-key map "\C-c\C-h" 'cocoa5-word-help)
    (define-key map "\C-c\C-e" 'cocoa5-close-block)
    (define-key map [f1] 'cocoa5-word-help) 
    (setq cocoa5-comint-mode-map map)))


;; ;; other useful things ?
;;   ( "\C-c\C-c" 'cocoa5-interrupt-or-comment-region)
;;   ( "\C-c\C-u" 'uncomment-region)
;;   ( "\C-c\C-e" 'cocoa5-close-block)
;;   ( "\C-c ]"   'cocoa5-close-block)
;;   ( "\C-c\C-o" 'cocoa5-toggle-selective-display-col)
;;   ( "\C-c\C-h" 'cocoa5-word-help)
;;   ( "\C-c\C-q" 'cocoa5-shell-quit)
;;   ( "\C-c\C-k" 'cocoa5-restart-shell)
;;   ( "\C-c\C-s" 'cocoa5-side-by-side)
;;   ( "\C-c\C-t" 'cocoa5-top-and-bottom)
;; )

;;;; fix: JAA 2016-06
(defun cocoa5-get-old-input ()
  "Corrected version of `comint-get-old-input-default'."
(interactive)
  (let (bof)
    (if (and (not comint-use-prompt-regexp)
             ;; Make sure we're in an input rather than output field.
             (null (get-char-property (setq bof (field-beginning)) 'field)))
	(field-string-no-properties bof)
      (comint-bol)
      (buffer-substring-no-properties (point)
				      (line-end-position)))))


;; ; jump to the error line.

;; ; ERROR: parse error in line 5 of device e:\burki\systems\cocoa\test.coc
;; ; -> put the cursor there.
;; ; ERROR: parse error in line 31 of device stdin
;; ; -> take no action.


(defun cocoa5-pop-to-buffer-source-find-X (err-or-warn)
  (progn 
    (setq found t)
    (setq etype err-or-warn)
    (setq emesg-start (point))
    (re-search-forward "^--> WHERE: " nil t)
    (beginning-of-line)
    (backward-char)
    (setq emesg (buffer-substring emesg-start (point)))
;;;;    (re-search-forward "^WHERE: at line(s?) \\([0-9]+\\) (column \\([0-9]+\\)) of \\(.*\\)$" nil t)
    (re-search-forward "^--> WHERE: at line\\(s?\\) \\([0-9]+\\) (column \\([0-9]+\\))" nil t)
    (progn
      (setq eline (buffer-substring (match-beginning 2) (match-end 2)))
      (setq ecol  (buffer-substring (match-beginning 3) (match-end 3)))
      (message "eline=%s, ecol=%s" eline ecol)
   ;  (setq efile (buffer-substring (match-beginning 3) (match-end 3)))
      (re-search-forward "of \\(.*\\)$" nil t)
      (setq erelfile (buffer-substring (match-beginning 1) (match-end 1)))
      (message "erelfile=%s" erelfile)
      (message "erelfile=%s, eline=%s, emesg=%s" erelfile eline emesg)
      )
    ;; (new) for adding path to "erelfile"
    (re-search-backward "source \"\\(.*\\)\";$" nil t)
    (setq efile (buffer-substring (match-beginning 1) (match-end 1)))
    )
)

(defun cocoa5-pop-to-buffer-source-find-error ()
  "If there is an error or warning, set the cursor at the error line."
  (interactive)
  (let ((found nil)
        (emesg nil) 
        (eline nil)
        (efile nil))
  (progn
;;    (message "debug cocoa5-find-error.")
;;    (switch-to-buffer "*cocoa5*") ; useless, since we are already there. 
    (cocoa5-pop-to-buffer-comint)
    (goto-char comint-last-input-start)
    ; search for an error message.
;    (if (re-search-forward "^Error \\(.*\\) at line \\([0-9]+\\) (column \\([0-9]+\\)) of \\(.*\\)$" nil t)
    (if (re-search-forward "--> ERROR: " nil t)
	(cocoa5-pop-to-buffer-source-find-X "ERROR")
      (if (re-search-forward "--> WARNING: " nil t)
	  (cocoa5-pop-to-buffer-source-find-X "WARNING")
	))
    (goto-char (point-max)) ; go to the end of the *cocoa5* buffer. 
    (message "DEBUG: cocoa5-find-error: found=%s, efile=%s, eline=%s, emesg=%s" found efile eline emesg)
;; AB22Apr2010: should not happen in cocoa-5
    ; messages such as "ERROR: parse error in line 31 of device stdin" 
    ; cannot be used for jumping to the error line.

    ; sometimes the file cannot be found (no abs path)
    ; in such a case, proceed as if no error was found.
    (if (and found (not (file-exists-p efile)))
      (progn
	(message "Emacs cannot find file \"%s\" (line %s col %s)\n%s: %s" efile eline ecol etype emesg)
	(setq found nil)
      ))
    (if found
      (progn 
         (cocoa5-pop-to-buffer-source) ;; <--- bad. ignores the efile. 
	 (find-file efile)            ;; <--- bad: doesn't behave like pop-to-buffer.

;;	 (goto-line (string-to-number eline))  AB22Apr2010: could have &optional BUFFER
	 (goto-line (string-to-number eline))
	 (forward-char (- (string-to-number ecol) 1))
	 (message "CoCoA %s (line %s col %s):\n%s" etype eline ecol emesg)
      )
      (cocoa5-pop-to-buffer-source) ; no error found. still, go back to the source buffer.
      )
   ))) 

(defun cocoa5-pop-to-buffer-source ()
  "Switch to the cocoa source buffer, if it exists."
  (interactive)
  (try-pop-to-buffer cocoa5-source-buffer))

(defun cocoa5-pop-to-buffer-comint ()
  "Switch to the *cocoa5* comint buffer, if it exists."
  (interactive)
  (try-pop-to-buffer "*cocoa5*"))

(defun try-pop-to-buffer (buffer)
  "Switch to buffer, if it exists."
  (let ((buf (get-buffer buffer)))
    (if buf
; switches buffer, but keeps the window layout:
; never changes the number and the sizes of the visible windows in the current frame.
	(pop-to-buffer buf)  
        ;; would destroy the window layout: (switch-to-buffer buf)
      (message "No buffer \"%s\"." buffer))))


; Optional feature: tab completion.
;; ANNA: achieved in cocoa5-electric-tab
;(add-hook 'cocoa5-mode-hook 
;   (lambda () (local-set-key (kbd "<tab>") 'dabbrev-expand)))

;; (add-hook 'cocoa5-comint-mode-hook 
;;    (lambda () (local-set-key (kbd "<tab>") 'dabbrev-expand)))

(defun cocoa5-syntax () 
  (progn 
  ;; Copied from the function cocoa5-mode. 
  ;; to do: avoid to duplicate code. 
  (setq local-abbrev-table cocoa5-mode-abbrev-table)
  (set-syntax-table cocoa5-mode-syntax-table)
  (make-local-variable 'indent-line-function)
  (setq indent-line-function 'cocoa5-indent-line)
  (make-local-variable 'case-fold-search)
  (setq case-fold-search t)

  (make-local-variable 'comment-start)
  (setq comment-start "--")
  (make-local-variable 'comment-end)
  (setq comment-end "")
  (make-local-variable 'comment-column)
  (setq comment-column 32)
  (make-local-variable 'comment-start-skip)
  (setq comment-start-skip  "/\\*+ *\\|// *")
  (make-local-variable 'comment-indent-function)
  (setq comment-indent-function 'cocoa5-indent-comment)
  (make-local-variable 'parse-sexp-ignore-comments)
  (setq parse-sexp-ignore-comments nil)

;; Font lock support
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults '(cocoa5-font-lock-keywords nil nil))
;; Imenu support
  (make-local-variable 'imenu-generic-expression)
  (setq imenu-generic-expression cocoa5-imenu-generic-expression)
;; Load wordlist
  (if cocoa5-auto-load-wordlist
	 (find-file-noselect cocoa5-wordlist-default-file 'NOWARN))
;; Menu
  (cocoa5-mode-setup-menubar)
  ))

; Optional feature: syntax highlighting. doesn't work fully. 
(add-hook 'cocoa5-comint-mode-hook 'cocoa5-syntax);; loses bindings in the menu!
(add-hook 'cocoa5-mode-hook 'cocoa5-syntax)

;(set 'pop-up-windows nil)

(defun cocoa5-exec-this-and-stay-in-comint (arg)
  "Send arg to CoCoA5
"
;  (cocoa5-shell-check-process)
;;  (cocoa5-shell-check-window) ; <---
  (save-excursion
;    (switch-to-buffer-other-window "*cocoa5*")
    (cocoa5-pop-to-buffer-comint)
;    (if (not (eq (process-status "cocoa5") 'run)) (cocoa5-make-shell))
    (goto-char (point-max))
;    (insert "debug ")
    (move-beginning-of-line nil)  (if (not (eq (point) (point-max))) (kill-line))
    (insert arg)
    (comint-send-input)
;    (switch-to-buffer-other-window cocoa5-file-buffer-name)
;    (other-window 1)

    ;; scroll to the beginning of the help
    (goto-char comint-last-input-start)
    ; todo: scroll so that the line with the cursor is the top line. 
    (recenter 0)

;    (insert "aasdf")
    )
;  (skip-chars-forward "-*\n-*")
)

;; ;(global-set-key [f1] 'cocoa5-word-help)

;; (defun cocoa-toggle-selective-display ()
;;   (interactive)
;;   (set-selective-display (if selective-display nil 1)))

(defun cocoa5-toggle-selective-display-col ()
  (interactive "*")
  (set-selective-display
   (if selective-display nil (or (1+ (current-column)) 1))))


;; ;; don't split windows:
;; ; (set 'pop-up-windows nil)

;; ;; don't create new frames:
;; ; (set 'pop-up-frames nil)
