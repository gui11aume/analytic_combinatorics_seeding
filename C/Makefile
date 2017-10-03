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

FIGS= simulp.pdf simulpdel.pdf simulpins.pdf simulp-inexact.pdf \
      simulp-average.pdf



all: $(FIGS) $(SIMULATIONS)

sim-sub: sim-sub.c mt.o
	$(CC) $(CFLAGS) sim-sub.c mt.o -o sim-sub

sim-ins: sim-ins.c mt.o
	$(CC) $(CFLAGS) sim-ins.c mt.o -o sim-ins

sim-del: sim-del.c mt.o
	$(CC) $(CFLAGS) sim-del.c mt.o -o sim-del

sim-inexact-sub: sim-inexact-sub.c mt.o
	$(CC) $(CFLAGS) sim-inexact-sub.c mt.o -o sim-inexact-sub

sim-average-sub: sim-average-sub.c mt.o
	$(CC) $(CFLAGS) sim-average-sub.c mt.o -o sim-average-sub

sim-average-del: sim-average-del.c mt.o
	$(CC) $(CFLAGS) sim-average-del.c mt.o -o sim-average-del

sim-average-del2: sim-average-del2.c mt.o
	$(CC) $(CFLAGS) sim-average-del2.c mt.o -o sim-average-del2

sim-fp: sim-fp.c mt.o
	$(CC) $(CFLAGS) sim-fp.c mt.o -o sim-fp

sim-MEM: sim-MEM.c mt.o
	$(CC) $(CFLAGS) sim-MEM.c mt.o -o sim-MEM

sim-MEM_fp: sim-MEM_fp.c mt.o
	$(CC) $(CFLAGS) sim-MEM_fp.c mt.o -o sim-MEM_fp



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
out-.08-average-sub.txt: sim-average-sub
	./sim-average-sub 50 17 0.08 > out-.08-average-sub.txt
	for ((i=52 ; i <= 150; i += 2)) ; do \
		./sim-average-sub $$i 17 0.08 \
			>> out-.08-average-sub.txt ; \
	done
out-.10-average-sub.txt: sim-average-sub
	./sim-average-sub 50 17 0.10 > out-.10-average-sub.txt
	for ((i=52 ; i <= 150; i += 2)) ; do \
		./sim-average-sub $$i 17 0.10 \
			>> out-.10-average-sub.txt ; \
	done
out-.12-average-sub.txt: sim-average-sub
	./sim-average-sub 50 17 0.12 > out-.12-average-sub.txt
	for ((i=52 ; i <= 150; i += 2)) ; do \
		./sim-average-sub $$i 17 0.12 \
			>> out-.12-average-sub.txt ; \
	done

# Average number of deletions
out-.05-.14-average-del.txt: sim-average-del
	./sim-average-del 50 17 0.05 0.14 > out-.05-.14-average-del.txt
	for ((i=52 ; i <= 150; i += 2)) ; do \
		./sim-average-del $$i 17 0.05 0.14 \
			>> out-.05-.14-average-del.txt ; \
	done
out-.05-.15-average-del.txt: sim-average-del
	./sim-average-del 50 17 0.05 0.15 > out-.05-.15-average-del.txt
	for ((i=52 ; i <= 150; i += 2)) ; do \
		./sim-average-del $$i 17 0.05 0.15 \
			>> out-.05-.15-average-del.txt ; \
	done
out-.05-.16-average-del.txt: sim-average-del
	./sim-average-del 50 17 0.05 0.16 > out-.05-.16-average-del.txt
	for ((i=52 ; i <= 150; i += 2)) ; do \
		./sim-average-del $$i 17 0.05 0.16 \
			>> out-.05-.16-average-del.txt ; \
	done

# Average number of deletions (2)
out-.05-.14-average-del2.txt: sim-average-del2
	./sim-average-del2 50 17 0.05 0.14 > out-.05-.14-average-del2.txt
	for ((i=52 ; i <= 150; i += 2)) ; do \
		./sim-average-del2 $$i 17 0.05 0.14 \
			>> out-.05-.14-average-del2.txt ; \
	done
out-.05-.15-average-del2.txt: sim-average-del2
	./sim-average-del2 50 17 0.05 0.15 > out-.05-.15-average-del2.txt
	for ((i=52 ; i <= 150; i += 2)) ; do \
		./sim-average-del2 $$i 17 0.05 0.15 \
			>> out-.05-.15-average-del2.txt ; \
	done
out-.05-.16-average-del2.txt: sim-average-del2
	./sim-average-del2 50 17 0.05 0.16 > out-.05-.16-average-del2.txt
	for ((i=52 ; i <= 150; i += 2)) ; do \
		./sim-average-del2 $$i 17 0.05 0.16 \
			>> out-.05-.16-average-del2.txt ; \
	done


# Substitutions (short)
out-.005.txt: sim-sub
	./sim-sub 25 17 0.005 > out-.005.txt
	for ((i=26 ; i <= 55; i += 1)) ; do \
		./sim-sub $$i 17 0.005 >> out-.005.txt ; \
	done
out-.010.txt: sim-sub
	./sim-sub 25 17 0.010 > out-.010.txt
	for ((i=26 ; i <= 55; i += 1)) ; do \
		./sim-sub $$i 17 0.010 >> out-.010.txt ; \
	done
out-.015.txt: sim-sub
	./sim-sub 25 17 0.015 > out-.015.txt
	for ((i=26 ; i <= 55; i += 1)) ; do \
		./sim-sub $$i 17 0.015 >> out-.015.txt ; \
	done

# False positives
out-.05-fp.txt: sim-fp
	./sim-fp 25 17 0.05 0.05 10000000 > out-.05-fp.txt
	for ((i=29 ; i <= 101; i += 4)) ; do \
		./sim-fp $$i 17 0.05 0.05 10000000 >> out-.05-fp.txt ; \
	done
	for ((i=103 ; i <= 151; i += 2)) ; do \
		./sim-fp $$i 17 0.05 0.05 100000000 >> out-.05-fp.txt ; \
	done
out-.15-fp.txt: sim-fp
	./sim-fp 25 17 0.05 0.15 10000000 > out-.15-fp.txt
	for ((i=29 ; i <= 101; i += 4)) ; do \
		./sim-fp $$i 17 0.05 0.15 10000000 >> out-.15-fp.txt ; \
	done
	for ((i=103 ; i <= 151; i += 2)) ; do \
		./sim-fp $$i 17 0.05 0.15 100000000 >> out-.15-fp.txt ; \
	done
out-.25-fp.txt: sim-fp
	./sim-fp 25 17 0.05 0.25 10000000 > out-.25-fp.txt
	for ((i=29 ; i <= 101; i += 4)) ; do \
		./sim-fp $$i 17 0.05 0.25 10000000 >> out-.25-fp.txt ; \
	done
	for ((i=103 ; i <= 151; i += 2)) ; do \
		./sim-fp $$i 17 0.05 0.25 100000000 >> out-.25-fp.txt ; \
	done


# MEMs
out-.08-MEM.txt: sim-MEM
	./sim-MEM 50 17 0.08 0.10 10000000 > out-.08-MEM.txt
	for ((i=52 ; i <= 200; i += 2)) ; do \
		./sim-MEM $$i 17 0.08 0.10 10000000 >> out-.08-MEM.txt ; \
	done
out-.10-MEM.txt: sim-MEM
	./sim-MEM 50 17 0.10 0.10 10000000 > out-.10-MEM.txt
	for ((i=52 ; i <= 200; i += 2)) ; do \
		./sim-MEM $$i 17 0.10 0.10 10000000 >> out-.10-MEM.txt ; \
	done
out-.12-MEM.txt: sim-MEM
	./sim-MEM 50 17 0.12 0.10 10000000 > out-.12-MEM.txt
	for ((i=52 ; i <= 200; i += 2)) ; do \
		./sim-MEM $$i 17 0.12 0.10 10000000 >> out-.12-MEM.txt ; \
	done

# False positives MEM
out-.05-MEM_fp.txt: sim-MEM_fp
	./sim-MEM_fp 25 17 0.05 0.05 10000000 > out-.05-MEM_fp.txt
	for ((i=29 ; i <= 101; i += 4)) ; do \
		./sim-MEM_fp $$i 17 0.05 0.05 10000000 \
			>> out-.05-MEM_fp.txt ; \
	done
	for ((i=103 ; i <= 151; i += 2)) ; do \
		./sim-MEM_fp $$i 17 0.05 0.05 100000000 \
			>> out-.05-MEM_fp.txt ; \
	done
out-.15-MEM_fp.txt: sim-MEM_fp
	./sim-MEM_fp 25 17 0.05 0.15 10000000 > out-.15-MEM_fp.txt
	for ((i=29 ; i <= 101; i += 4)) ; do \
		./sim-MEM_fp $$i 17 0.05 0.15 10000000 \
			>> out-.15-MEM_fp.txt ; \
	done
	for ((i=103 ; i <= 151; i += 2)) ; do \
		./sim-MEM_fp $$i 17 0.05 0.15 100000000 \
			>> out-.15-MEM_fp.txt ; \
	done
out-.25-MEM_fp.txt: sim-MEM_fp
	./sim-MEM_fp 25 17 0.05 0.25 10000000 > out-.25-MEM_fp.txt
	for ((i=29 ; i <= 101; i += 4)) ; do \
		./sim-MEM_fp $$i 17 0.05 0.25 10000000 \
			>> out-.25-MEM_fp.txt ; \
	done
	for ((i=103 ; i <= 151; i += 2)) ; do \
		./sim-MEM_fp $$i 17 0.05 0.25 100000000 \
			>> out-.25-MEM_fp.txt ; \
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
	R -f make_figure_simulation_average_sub.R

simuldel-average.pdf: out-.05-.14-average-del.txt \
		out-.05-.15-average-del.txt out-.05-.16-average-del.txt
	R -f make_figure_simulation_average_del.R

simuldel-average2.pdf: out-.05-.14-average-del2.txt \
		out-.05-.15-average-del2.txt out-.05-.16-average-del2.txt
	R -f make_figure_simulation_average_del2.R

simulp_short.pdf: out-.005.txt out-.010.txt out-.015.txt
	R -f make_figure_simulation_sub_short.R

simulp_approx.pdf: out-.05-.15-.04-.45.txt out-.05-.15-.05-.45.txt \
		out-.05-.15-.06-.45.txt
	R -f make_figure_simulation_sub_approx.R

simulp_false_positives.pdf: out-.05-fp.txt out-.15-fp.txt out-.25-fp.txt
	R -f make_figure_simulation_false_positives.R

simulp_mem.pdf: out-.08-MEM.txt out-.10-MEM.txt out-.12-MEM.txt
	R -f make_figure_simulation_mem.R

simulp_mem_fp.pdf: out-.05-MEM_fp.txt out-.15-MEM_fp.txt out-.25-MEM_fp.txt
	R -f make_figure_simulation_mem_fp.R


clean:
	rm -f $(SIMULATIONS) $(EXEC) $(FIGS)
