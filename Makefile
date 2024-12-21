CFILES=main.c
CFLAGS=\
	-O3 -lm -lraylib -fsanitize=address \
	-I/opt/homebrew/Cellar/raylib/5.5/include -L/opt/homebrew/Cellar/raylib/5.5/lib

same: $(CFILES)
	gcc main.c -o same $(CFLAGS)

run: same
	./same

clean:
	rm -f same

.PHONY: run clean