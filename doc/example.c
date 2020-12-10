
#include <swarm.h>
#include <swarm_random.h>

static void routine (THREADED)
{
	int i, n;
	int elems_per_thread = 5, max_random = 10, partial_sum = 0, prefix_sum = 0, global_sum = 0;
	int *arr;
	double global_mean = 0, partial_squared_sum = 0, squared_sum = 0, std_dev = 0;

	n = elems_per_thread * THREADS;

	arr = SWARM_malloc(n * sizeof(int), TH);
	on_one_thread
	{
		printf("Array memory allocated on a single thread...\n\n");
	}

	SWARM_random_init(TH);
    SWARM_srandom(MYTHREAD + 1,TH);
	
	pardo(i, 0, n, 1)
	{
		arr[i] = SWARM_random(TH)%max_random;
	}

	SWARM_Barrier();

	on_one_thread
	{
		printf("Randomly generated array is:\n");
		for(i = 0; i < n; i++)
		{
			printf("arr[%d] = %ld\n", i, arr[i]);
		}
		printf("\n\n");
	}

	SWARM_Barrier();
	
	for(i = 0; i < elems_per_thread; i++)
	{
		partial_sum += arr[i + MYTHREAD*elems_per_thread];
	}

	on_one_thread
	{
		printf("Partial Sum:\n");
	}
	
	printf("Thread %d: %d\n", MYTHREAD, partial_sum);

	SWARM_Barrier();

	on_one_thread
	{
		printf("\n\n");
		printf("Global sum calculated using Reduce operation...\n\n");
		printf("Global Sum:\n");
	}
	
	SWARM_Barrier();
	global_sum = SWARM_Reduce_i(partial_sum, SUM, TH);
	printf("Thread %d: %d\n", MYTHREAD, global_sum);

	on_one_thread
	{
		global_mean = 1.0 * global_sum/n;
		printf("\n\n");
		printf("Mean calculated on a thread = %f\n", global_mean);
		printf("\n\n");
		printf("Broadcasting mean...\n\n");
		printf("Global Mean:\n");
	}

	global_mean = SWARM_Bcast_d(global_mean, TH);

	SWARM_Barrier();
	printf("Thread %d: %f\n", MYTHREAD, global_mean);
	SWARM_Barrier();

	for(i = 0; i < elems_per_thread; i++)
	{
		partial_squared_sum += (arr[i + MYTHREAD*elems_per_thread] - global_mean)*(arr[i + MYTHREAD*elems_per_thread] - global_mean);
	}

	SWARM_Barrier();
	
	on_one_thread
	{
		printf("\n\n");
		printf("Partial Squared Sum:\n");
	}
	SWARM_Barrier();
	printf("Thread %d: %f\n", MYTHREAD, partial_squared_sum);
	SWARM_Barrier();

	on_one_thread
	{
		printf("\n\nCalculating total squared sum from partial squared sums using Scan operation...\n\n");
		printf("Total Squared Sum:\n");
	}

	SWARM_Barrier();
	squared_sum = SWARM_Scan_d(partial_squared_sum, SUM, TH);
		
	printf("Thread %d: %f\n", MYTHREAD, squared_sum);
	SWARM_Barrier();

	on_thread(THREADS - 1)
	{
		printf("\n\n");
		printf("Calculating Standard Deviation on last thread...\n\n");
		std_dev = sqrt(squared_sum / n);
		printf("Standard Deviation = %f\n", std_dev);
	}

	SWARM_Barrier();
	
	on_one_thread
	{
		printf("\n\n");
		printf("Releasing array memory...\n\n");
	}
	SWARM_free(arr, TH);
}

int main (int argc, char **argv)
{
     SWARM_Init(&argc,&argv);

     SWARM_Run ((void *)routine);

     SWARM_Finalize();

     return 0;
}
