SHELL= bash
go: sim
	./sim 100 17 0.05 0.14 > out-.05-.14.txt
	for ((i=104 ; i <= 350; i += 4)) ; do \
		./sim $$i 17 0.05 0.14 >> out-.05-.14.txt ; \
	done
	for ((i=352 ; i <= 500; i += 2)) ; do \
		./sim $$i 17 0.05 0.14 >> out-.05-.14.txt ; \
	done
#	./sim 100 17 0.05 0.15 > out-.05-.15.txt
#	for ((i=104 ; i <= 350; i += 4)) ; do \
#		./sim $$i 17 0.05 0.15 >> out-.05-.15.txt ; \
#	done
#	for ((i=352 ; i <= 500; i += 2)) ; do \
#		./sim $$i 17 0.05 0.15 >> out-.05-.15.txt ; \
#	done
	./sim 100 17 0.05 0.16 > out-.05-.16.txt
	for ((i=104 ; i <= 350; i += 4)) ; do \
		./sim $$i 17 0.05 0.16 >> out-.05-.16.txt ; \
	done
	for ((i=352 ; i <= 500; i += 2)) ; do \
		./sim $$i 17 0.05 0.16 >> out-.05-.16.txt ; \
	done

sim: simdel.c
	gcc -std=gnu99 -g -O3 simdel.c -o sim

clean:
	rm go
