export ROOT = $(realpath .)
export OBJ_DIR = ${ROOT}/obj/

include inc.mk

LIST = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50

BIN_DIR = ./bin/
TARGET := $(BIN_DIR)/main

SUB_DIRS := prf 
SUB_DIRS += decode 
SUB_DIRS += hash 
SUB_DIRS += gf2x
SUB_DIRS += common

CSRC = kem.c

OBJS = $(OBJ_DIR)/*.o
ifdef USE_NIST_RAND
    CSRC += FromNIST/rng.c FromNIST/PQCgenKAT_kem.c
    OBJS += $(OBJ_DIR)/FromNIST/*.o
else
    SUB_DIRS += tests
endif

.PHONY: $(SUB_DIRS)

include rules.mk

SRC_FOR_TIDY := ${ROOT}/*.c
SRC_FOR_TIDY += prf/*.c
SRC_FOR_TIDY += decode/*.c
SRC_FOR_TIDY += hash/*.c
SRC_FOR_TIDY += gf2x/*.c
SRC_FOR_TIDY += common/*.c

all: $(BIN_DIR) $(OBJ_DIR) $(SUB_DIRS)
	$(CC) $(OBJS) $(CFLAGS) $(EXTERNAL_LIBS) -o $(TARGET) -lm

$(SUB_DIRS):
	make -C $@

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)
	mkdir -p $(OBJ_DIR)/FromNIST

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -rf $(OBJ_DIR)
	rm -rf $(BIN_DIR)

pretty:
	find . -name '*.c' -o -name '*.h' | xargs clang-format-9 -style=file -i

tidy:
	clang-tidy-9 ${SRC_FOR_TIDY} -p $(ROOT) --fix-errors --format-style=file -- ${CFLAGS}

run:
	cd bin;\
	for i in $(LIST); do \
		mkdir $$i;\
		cp main $$i;\
		cd $$i;\
		nohup ./main &\
		cd ..;\
		sleep 1;\
		echo "第 $$i 个进程生成中...";\
	done
	@echo "---运行完成, 所有程序将在后台执行---"

copy:
	cd bin;\
	mkdir all;\
	for i in $(LIST); do \
		cd $$i;\
		cp iter_data_all.txt ../all;\
		cd ..;\
		cd all;\
		mv iter_data_all.txt iter_data_all_$$i;\
		cd ..;\
	done
	cp copy_data/copy.c bin/all;\
	cd bin/all;\
	gcc copy.c -o copy;\
	./copy;\
	cd ..;\
	mkdir iter_data_all;\
	cd all;\
	cp iter_data_all ../iter_data_all;\
	cd ..;\
	cd iter_data_all;\
	wc iter_data_all

cleans:
	cd bin;\
	rm -rf iter_data_all all 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64
