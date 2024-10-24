FC = ifort
OBJS = main.o calc_force.o velocity.o verlet.o initpos.o \
       gaus_setup.o gaus_update.o add_gaus.o velocity_mod.o neb.o\

PROGRAM = md

.f.o:
	$(FC) -c $< -o $@
.F.o:
	$(FC) -FR -traceback -fast -tpp7 -c $< -o $@

# Default
$(PROGRAM) : $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	$(FC) $(CFLAGS) $(OBJS) -o $(PROGRAM)
configmake : configmake.o
	@echo "Linking $(PROGRAM) ..."
	$(FC) $(CFLAGS) configmake.o -o configmake

clean :; 
	rm -rf $(OBJS)
	rm -rf $(PROGRAM)
	rm -rf configmake configmake.o
	rm -rf *.mod
	rm -rf *~
