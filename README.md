# HW3_APD

A Project ilustrating the usefulness of MPI in Distributed Programming.

## Implementation

Initially, the leaders obtain their worker processes from the input file to form a special vector of vectors, namely the topology. Obviously, each Cluster object retains its own workers in the private vector of the class.

Of course, the leaders send their rank so that the worker processes know how to communicate with them. Moreover, the leaders send their own worker processes to the right side and receive worker processes belonging to other leaders from the left side so that, in the end, the topology is known. Any process will display the topology when it knows it.

At the beginning, leader 0 generates the vector and sends it:
- to its own worker processes.
- to the right side.

This method is repeated among the leaders, the final goal being to obtain the calculated vector in process 0.

The worker processes obtain the vector and handle it based on the relationship in the document. They will pass the processed vector to their leader. Obviously, all intermediate processes will transfer the obtained vector to the left side (to reach proc. 0).

Towards the end, leader 0 displays the desired vector obtained with the help of the 'recv()' function from MPI.

The parts related to the study of topology and vector calculation are approximately similar, a difference being the way of traversing the graph / topology. More precisely, leader 0 can only send / receive from the right side, leaders 2 and 3 from both sides, and leader 1 only from the left side. Regarding that vector, the route changes as follows: from proc 0 -> proc 1 -> proc 2 -> proc 3 to proc 0 -> proc 3 -> proc 2 -> proc 1.

## Bonus

For the bonus, an attempt was made to isolate process 1 from the others. It can only communicate with its own worker processes.

Obviously, if leader 1 is ignored, then its children will also bear this. Given that 1 no longer participates, the rank of the current process and its number of workers are calculated in order to use the same relationship from the problem statement.

To determine the aforementioned, the topology studied by each participating process is traversed.

Some references that helped with debugging:
- https://education.molssi.org/parallel-programming/04-distributed-examples/index.html
- https://www.open-mpi.org/faq/?category=troubleshooting
- https://rantahar.github.io/introduction-to-mpi/aio/index.html
- https://cplusplus.com/forum/general/282234/
