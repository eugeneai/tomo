;--Build the Y array
    ;***
    IF VERBOSE THEN PRINTW, WID=verbose,  function_string+' Building array Y'
    starttime =  SYSTIME(/SECONDS)

    muforWterm = 0
    muforVterm = 0

    FOR n=1,Nmax  DO BEGIN
    FOR l=0,Lmax  DO BEGIN
    FOR m=0,l     DO BEGIN


       ;--Calculate the W terms (Real Cosine) and the V terms (Imaginary Sine)
       ;
       IF KEYWORD_SET(fast) THEN Term = YarrayTerm_FAST(p,alp,bet,gam, gbar, precalculated, param=[l,m,n]) $
       ELSE Term   = YarrayTerm(p,alp,bet,gam, gbar, param=[l,m,n])

       FOR i=1,Ndim     DO BEGIN

         ;--Place the W term and V terms at their respective positions in the Y array
         ;
                   muforWterm  = (Lmax+1)^2*( (i-1)*Nmax + n - 1 ) + l^2 + m + 1
         IF m NE 0 THEN  muforVterm  = (Lmax+1)^2*( (i-1)*Nmax + n - 1  ) + l^2 + l + m + 1

                   Yarray[muforWterm-1] = Term[0,i-1]
         IF m NE 0 THEN  Yarray[muforVterm-1] = Term[1,i-1]



       ENDFOR ;--end loop over i

    ENDFOR ;--end loop over m
    ENDFOR ;--end loop over l
    ENDFOR ;--end loop over n


