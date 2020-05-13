CC= g++
CFLAGS=	-O3 -Wall
FILES=	accept1 accept2 accept3 btraj grang2 grang grtheta \
	gr_uh hoh ouo tiltd tiltx traj trajedy



all:    $(FILES)
	$(CC) $(CFLAGS) $<.c -o $@