OBJ_FILES += $(patsubst %.S, $(OBJ_DIR)/%.o, $(SSRC))
OBJ_FILES += $(patsubst %.c, $(OBJ_DIR)/%.o, $(CSRC))

all: $(OBJ_FILES)

$(OBJ_DIR)/%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: FromNIST/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: %.S
	$(CC) $(CFLAGS) -c -o $@ $<
