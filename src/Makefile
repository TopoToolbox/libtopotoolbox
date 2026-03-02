.POSIX:
.SUFFIXES:

SRCS=hillshade.c drainagebasins.c knickpoints.c excesstopography.c fillsinks.c flow_accumulation.c flow_routing.c gradient8.c gwdt.c identifyflats.c reconstruct.c streamquad.c topotoolbox.c swaths.c graphflood/gf_utils.c graphflood/sfgraph.c graphflood/priority_flood_standalone.c graphflood/gf_flowacc.c graphflood/graphflood.c helpers/priority_queue.c helpers/dijkstra.c helpers/polyline.c helpers/stat_func.c helpers/deque.c

OBJS=$(SRCS:.c=.o)

libtopotoolbox.a: $(OBJS)
	$(AR) $(ARFLAGS) libtopotoolbox.a $(OBJS)
	rm -f $(OBJS)

.SUFFIXES: .c .o
.c.o:
	$(CC) $(CFLAGS) -DTOPOTOOLBOX_STATICLIB -DNDEBUG -I../include -c -o $@ $<