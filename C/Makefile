SHELL= bash
CC= gcc
CFLAGS= -std=gnu99 -O3

EXEC= sim-sub sim-ins sim-del

SIMULATIONS=                       \
	out-.08.txt                \
	out-.10.txt                \
	out-.12.txt                \
	out-.05-.14.txt            \
	out-.05-.15.txt            \
	out-.05-.16.txt            \
	out-.05-.15-.04-.45.txt    \
	out-.05-.15-.05-.45.txt    \
	out-.05-.15-.06-.45.txt    \
	out-.08-inexact.txt        \
	out-.10-inexact.txt        \
	out-.12-inexact.txt        \
	out-.08-average.txt        \
	out-.10-average.txt        \
	out-.12-average.txt

FIGS= simulp.pdf simulpdel.pdf simulpins.pdf simulp-inexact.pdf



all: $(FIGS) $(SIMULATIONS)

sim-sub: sim-sub.c mt.o
	$(CC) $(CFLAGS) sim-sub.c mt.o -o sim-sub

sim-ins: sim-ins.c mt.o
	$(CC) $(CFLAGS) sim-ins.c mt.o -o sim-ins

sim-del: sim-del.c mt.o
	$(CC) $(CFLAGS) sim-del.c mt.o -o sim-del

sim-inexact-sub: sim-inexact-sub.c mt.o
	$(CC) $(CFLAGS) sim-inexact-sub.c mt.o -o sim-inexact-sub

sim-average: sim-average.c mt.o
	$(CC) $(CFLAGS) sim-average.c mt.o -o sim-average



# Substitutions.
out-.08.txt: sim-sub
	./sim-sub 50 17 0.08 > out-.08.txt
	for ((i=52 ; i <= 200; i += 2)) ; do \
		./sim-sub $$i 17 0.08 >> out-.08.txt ; \
	done
out-.10.txt: sim-sub
	./sim-sub 50 17 0.10 > out-.10.txt
	for ((i=52 ; i <= 200; i += 2)) ; do \
		./sim-sub $$i 17 0.10 >> out-.10.txt ; \
	done
out-.12.txt: sim-sub
	./sim-sub 50 17 0.12 > out-.12.txt
	for ((i=52 ; i <= 200; i += 2)) ; do \
		./sim-sub $$i 17 0.12 >> out-.12.txt ; \
	done

# Substitutions and deletions.
out-.05-.14.txt: sim-del
	./sim-del 100 17 0.05 0.14 > out-.05-.14.txt
	for ((i=104 ; i <= 350; i += 4)) ; do \
		./sim-del $$i 17 0.05 0.14 >> out-.05-.14.txt ; \
	done
	for ((i=352 ; i <= 500; i += 2)) ; do \
		./sim-del $$i 17 0.05 0.14 >> out-.05-.14.txt ; \
	done
out-.05-.15.txt: sim-del
	./sim-del 100 17 0.05 0.15 > out-.05-.15.txt
	for ((i=104 ; i <= 350; i += 4)) ; do \
		./sim-del $$i 17 0.05 0.15 >> out-.05-.15.txt ; \
	done
	for ((i=352 ; i <= 500; i += 2)) ; do \
		./sim-del $$i 17 0.05 0.15 >> out-.05-.15.txt ; \
	done
out-.05-.16.txt: sim-del
	./sim-del 100 17 0.05 0.16 > out-.05-.16.txt
	for ((i=104 ; i <= 350; i += 4)) ; do \
		./sim-del $$i 17 0.05 0.16 >> out-.05-.16.txt ; \
	done
	for ((i=352 ; i <= 500; i += 2)) ; do \
		./sim-del $$i 17 0.05 0.16 >> out-.05-.16.txt ; \
	done

# Substitutions deletions and insertions.
out-.05-.15-.04-.45.txt: sim-ins
	./sim-ins 100 17 0.05 0.15 0.04 0.45 > out-.05-.15-.04-.45.txt
	for ((i=108 ; i <= 604; i += 12)) ; do \
	   ./sim-ins $$i 17 0.05 0.15 0.04 0.45 >> out-.05-.15-.04-.45.txt ; \
	done
	for ((i=612 ; i <= 1000; i += 8)) ; do \
	   ./sim-ins $$i 17 0.05 0.15 0.04 0.45 >> out-.05-.15-.04-.45.txt ; \
	done
out-.05-.15-.05-.45.txt: sim-ins
	./sim-ins 100 17 0.05 0.15 0.05 0.45 > out-.05-.15-.05-.45.txt
	for ((i=108 ; i <= 604; i += 12)) ; do \
	   ./sim-ins $$i 17 0.05 0.15 0.05 0.45 >> out-.05-.15-.05-.45.txt ; \
	done
	for ((i=612 ; i <= 1000; i += 8)) ; do \
	   ./sim-ins $$i 17 0.05 0.15 0.05 0.45 >> out-.05-.15-.05-.45.txt ; \
	done
out-.05-.15-.06-.45.txt: sim-ins
	./sim-ins 100 17 0.05 0.15 0.06 0.45 > out-.05-.15-.06-.45.txt
	for ((i=108 ; i <= 604; i += 12)) ; do \
	   ./sim-ins $$i 17 0.05 0.15 0.06 0.45 >> out-.05-.15-.06-.45.txt ; \
	done
	for ((i=612 ; i <= 1000; i += 8)) ; do \
	   ./sim-ins $$i 17 0.05 0.15 0.06 0.45 >> out-.05-.15-.06-.45.txt ; \
	done

# Substitutions (inexact)
out-.08-inexact.txt: sim-inexact-sub
	./sim-inexact-sub 50 17 0.08 10000000 > out-.08-inexact.txt
	for ((i=52 ; i <= 148; i += 2)) ; do \
		./sim-inexact-sub $$i 17 0.08 10000000 \
			>> out-.08-inexact.txt ; \
	done
	for ((i=150 ; i <= 200; i += 2)) ; do \
		./sim-inexact-sub $$i 17 0.08 100000000 \
			>> out-.08-inexact.txt ; \
	done
out-.10-inexact.txt: sim-inexact-sub
	./sim-inexact-sub 50 17 0.10 10000000 > out-.10-inexact.txt
	for ((i=52 ; i <= 148; i += 2)) ; do \
		./sim-inexact-sub $$i 17 0.10 10000000 \
			>> out-.10-inexact.txt ; \
	done
	for ((i=150 ; i <= 200; i += 2)) ; do \
		./sim-inexact-sub $$i 17 0.10 100000000 \
			>> out-.10-inexact.txt ; \
	done
out-.12-inexact.txt: sim-inexact-sub
	./sim-inexact-sub 50 17 0.12 10000000 > out-.12-inexact.txt
	for ((i=52 ; i <= 148; i += 2)) ; do \
		./sim-inexact-sub $$i 17 0.12 10000000 \
		>> out-.12-inexact.txt ; \
	done
	for ((i=150 ; i <= 200; i += 2)) ; do \
		./sim-inexact-sub $$i 17 0.12 100000000 \
		>> out-.12-inexact.txt ; \
	done

# Average number of substitutions
out-.08-average.txt: sim-average
	./sim-average 50 17 0.08 > out-.08-average.txt
	for ((i=52 ; i <= 150; i += 2)) ; do \
		./sim-average $$i 17 0.08 \
			>> out-.08-average.txt ; \
	done
out-.10-average.txt: sim-average
	./sim-average 50 17 0.10 > out-.10-average.txt
	for ((i=52 ; i <= 150; i += 2)) ; do \
		./sim-average $$i 17 0.10 \
			>> out-.10-average.txt ; \
	done
out-.12-average.txt: sim-average
	./sim-average 50 17 0.12 > out-.12-average.txt
	for ((i=52 ; i <= 150; i += 2)) ; do \
		./sim-average $$i 17 0.12 \
			>> out-.12-average.txt ; \
	done


simulp.pdf: out-.08.txt out-.10.txt out-.12.txt
	R -f make_figure_simulation_sub.R

simulpdel.pdf: out-.05-.14.txt out-.05-.15.txt out-.05-.16.txt
	R -f make_figure_simulation_del.R

simulpins.pdf: out-.05-.15-.04-.45.txt out-.05-.15-.05-.45.txt \
		out-.05-.15-.06-.45.txt
	R -f make_figure_simulation_ins.R

simulp-inexact.pdf: out-.08-inexact.txt out-.10-inexact.txt \
		out-.12-inexact.txt
	R -f make_figure_simulation_inexact_sub.R

simulp-average.pdf: out-.08-average.txt out-.10-average.txt \
		out-.12-average.txt
	R -f make_figure_simulation_average.R


clean:
	rm -f $(SIMULATIONS) $(EXEC) $(FIGS)
