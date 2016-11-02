

# original source files:
#SRCS = disnan.f dlabad.f dladiv.f dlaisnan.f dlamch.f dlapy3.f dznrm2.f ieeeck.f ilaenv.f ilazlc.f ilazlr.f iparmq.f lsame.f xerbla.f zcopy.f zdscal.f zgelq2.f zgelqf.f zgels.f zgemm.f zgemv.f zgeqr2.f zgeqrf.f zgerc.f zlacgv.f zladiv.f zlange.f zlarf.f zlarfb.f zlarfg.f zlarft.f zlascl.f zlaset.f zlassq.f zscal.f ztrmm.f ztrmv.f ztrsm.f ztrtrs.f zunm2r.f zunml2.f zunmlq.f zunmqr.f 

# original object files:
#OBJS = disnan.o dlabad.o dladiv.o dlaisnan.o dlamch.o dlapy3.o dznrm2.o ieeeck.o ilaenv.o ilazlc.o ilazlr.o iparmq.o lsame.o xerbla.o zcopy.o zdscal.o zgelq2.o zgelqf.o zgels.o zgemm.o zgemv.o zgeqr2.o zgeqrf.o zgerc.o zlacgv.o zladiv.o zlange.o zlarf.o zlarfb.o zlarfg.o zlarft.o zlascl.o zlaset.o zlassq.o zscal.o ztrmm.o ztrmv.o ztrsm.o ztrtrs.o zunm2r.o zunml2.o zunmlq.o zunmqr.o 

OBJS = disnanquad.o dlabadquad.o dladivquad.o dlaisnanquad.o dlamchquad.o dlapy3quad.o dznrm2quad.o ieeeckquad.o ilaenvquad.o ilazlcquad.o ilazlrquad.o iparmqquad.o lsamequad.o xerblaquad.o zcopyquad.o zdscalquad.o zgelq2quad.o zgelqfquad.o zgelsquad.o zgemmquad.o zgemvquad.o zgeqr2quad.o zgeqrfquad.o zgercquad.o zlacgvquad.o zladivquad.o zlangequad.o zlarfquad.o zlarfbquad.o zlarfgquad.o zlarftquad.o zlasclquad.o zlasetquad.o zlassqquad.o zscalquad.o ztrmmquad.o ztrmvquad.o ztrsmquad.o ztrtrsquad.o zunm2rquad.o zunml2quad.o zunmlqquad.o zunmqrquad.o 


F77 = mpif90
#F77 = gfortran
#F77 = ifort
F77FLAGS =  -O2

GENLIB= libzgelsquad.a

   
# $(PROG): $(OBJS)

#libzgelsquad.a: $(OBJS)
#	ar r libzgelsquad.a $(OBJS)
   
$(GENLIB): $(OBJS)
	ar r $(GENLIB) $(OBJS)

    #$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)


.SUFFIXES: $(SUFFIXES) .f

.f.o:
	$(F77) $(F77FLAGS) -c $<



clean:
	rm -f $(GENLIB) $(OBJS) *.mod

