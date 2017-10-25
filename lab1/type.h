typedef struct node {
	int id;
	struct node* next;
}node;
typedef struct box {
	int id;
	int location[4];
	int topnum_neigh;
    struct node* toplist;
	int bottomnum_neigh;
	struct node* bottomlist;
	int leftnum_neigh;
	struct node* leftlist;
	int rightnum_neigh;
	struct node* rightlist;
	float temperature;
	float adjacent_temp;
} box;
//int clock_gettime(clockid_t clk_id, struct timespec *tp);