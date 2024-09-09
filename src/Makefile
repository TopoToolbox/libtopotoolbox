.POSIX:
.SUFFIXES:

SRCS=excesstopography.c fillsinks.c flow_accumulation.c flow_routing.c gradient8.c gwdt.c identifyflats.c morphology/reconstruct.c priority_queue.c topotoolbox.c graphflood/gf_utils.c graphflood/sfgraph.c graphflood/priority_flood_standalone.c graphflood/gf_flowacc.c graphflood/graphflood.c

OBJS=$(SRCS:.c=.o)

libtopotoolbox.a: $(OBJS)
	$(AR) $(ARFLAGS) libtopotoolbox.a $(OBJS)
	rm -f $(OBJS)

.SUFFIXES: .c .o
.c.o:
	$(CC) $(CFLAGS) -DTOPOTOOLBOX_STATICLIB -I../include -c -o $@ $<
