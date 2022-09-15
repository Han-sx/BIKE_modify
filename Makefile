export ROOT = $(realpath .)
export OBJ_DIR = ${ROOT}/obj/

include inc.mk

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
	$(CC) $(OBJS) $(CFLAGS) $(EXTERNAL_LIBS) -o $(TARGET)

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

