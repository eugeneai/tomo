--Each term in the Y array
;  Balandin code: RIGHTPART_WWVEC, l.1288
;  (identical to Balandin's method)
;***
FUNCTION YarrayTerm, p, alp, bet, gam, gbar, param=param
    Nump    = N_ELEMENTS(p)
    Numalp     = N_ELEMENTS(alp)
    Numbet     = N_ELEMENTS(bet)
    Numgam     = N_ELEMENTS(gam)
    dp         = p[1] - p[0]
    dalp       = alp[1] - alp[0]
    dbet     = bet[1] - bet[0]
    dgam     = gam[1] - gam[0]
    l       = param[0] &   m  = param[1] &   n  = param[2]
    fullintegral= FLTARR(2, 3)  ; [ 0=W term 1=V term  ; i=1,2,3 ]


;   ;--Balandin's routine has a different array for p (see l.3069 vs l. 3083)
;   ;  (doesn't make a big difference)
;   RRE    = 1.0
;   DD     = 1.5
;   bkappa0   = ASIN(RRE/DD)
;   dkappa    = 2.*bkappa0/(Nump-1)
;   bkapparange = [-bkappa0, bkappa0]
;   bkappa    = FINDGEN(Nump)/(Nump-1)*(bkapparange[1]-bkapparange[0]) + bkapparange[0]
;   pp        = DD*SIN(bkappa)
;   p=pp

    ;==integral, over beta first,
    FOR ixb = 0, Numbet-1 DO BEGIN
       beta      = bet[ixb]
       Sinbeta   = SIN(beta)
       wb       = (ixb EQ 0 OR ixb EQ Numbet-1) ? 0.5 : 1.0        ;--trapezium rule integration
       gintegrand    = fullintegral*0.0
       ;==integral, over gamma,
       FOR ixg = 0, Numgam-1 DO BEGIN
         gammaa    = gam[ixg]
         wg         = (ixg EQ 0 OR ixg EQ Numgam-1) ? 0.5 : 1.0        ;--trapezium rule integration
         aintegrand   = fullintegral*0.0
         ;==integral, over alpha,
         FOR ixa = 0, Numalp-1 DO BEGIN
          alpha    = alp[ixa]
          wa         = (ixa EQ 0 OR ixa EQ Numalp-1) ? 0.5 : 1.0    ;--trapezium rule integration
          nv          = NVEC(alpha,beta,gammaa);REFORM(Euler3DRotMatrix(alpha,beta,gammaa) # [0.0, 0.0, 1.0] )
          pintegrand  = fullintegral*0.0
          ;==integral, over p range, of gbar*Wfunc
          FOR ixp = 0, Nump-1   DO BEGIN
                        pnt    = p[ixp]
                        wp        = (ixp EQ 0 OR ixp EQ Nump-1) ? 0.5 : 1.0        ;--trapezium rule integration
                                 pintegrand[0,*] += gbar[ixp,ixa,ixb,ixg] * Wfunc(pnt,alpha,beta,gammaa, param=[l,m,n]) * dp * wp
                 IF m NE 0 THEN    pintegrand[1,*] += gbar[ixp,ixa,ixb,ixg] * Vfunc(pnt,alpha,beta,gammaa, param=[l,m,n]) * dp * wp
              ENDFOR
              ;==end loop over p
          aintegrand[*,0] += pintegrand * nv[0] * dalp * wa
          aintegrand[*,1] += pintegrand * nv[1] * dalp * wa
          aintegrand[*,2] += pintegrand * nv[2] * dalp * wa
         ENDFOR
         ;==end loop over alpha
         gintegrand += aintegrand * dgam * wg
       ENDFOR
       ;==end loop over gamma
       fullintegral += gintegrand * Sinbeta * dbet * wb
    ENDFOR
    ;==end loop over beta

    RETURN, fullintegral
END