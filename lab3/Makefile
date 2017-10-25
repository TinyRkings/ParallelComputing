all: disposable persistent
.PHONY : clean 
disposable: lou_shuangsheng_disposable.o
	gcc lou_shuangsheng_disposable.o -o disposable -lrt -O3 -fopenmp
persistent: lou_shuangsheng_persistent.o
	gcc lou_shuangsheng_persistent.o -o persistent -lrt -O3 -fopenmp
lou_shuangsheng_disposable.o: lou_shuangsheng_disposable.c
	gcc -c lou_shuangsheng_disposable.c -lrt -O3 -fopenmp
lou_shuangsheng_persistent.o: lou_shuangsheng_persistent.c
	gcc -c lou_shuangsheng_persistent.c -lrt -O3 -fopenmp

clean:
	rm -rf *.o *.out disposable persistent