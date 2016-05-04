(custom-field-function/define
 '(((name gsd) (display "e ^ (sqrt udm-2)") (syntax-tree ("**" 2.71828 ("sqrt" "udm-2"))) (code (field-** 2.71828 (field-sqrt (field-load "udm-2")))))
   ((name ntot) (display "udm-0") (syntax-tree "udm-0") (code (field-load "udm-0")))
   ((name cmd-nm) (display "udm-1 * 10 ^ 9") (syntax-tree ("*" "udm-1" 1000000000.)) (code (field-* (field-load "udm-1") 1000000000.)))
   ))