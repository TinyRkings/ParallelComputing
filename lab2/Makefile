all: disposable persistent
.PHONY : clean 
disposable: lou_shuangsheng_disposable.o
	gcc lou_shuangsheng_disposable.o -o disposable -lpthread -lrt -O3
persistent: lou_shuangsheng_persistent.o
	gcc lou_shuangsheng_persistent.o -o persistent -lpthread -lrt -O3
lou_shuangsheng_disposable.o: lou_shuangsheng_disposable.c
	gcc -c lou_shuangsheng_disposable.c
lou_shuangsheng_persistent.o: lou_shuangsheng_persistent.c
	gcc -c lou_shuangsheng_persistent.c

clean:
	rm -rf *.o *.out disposable persistent