# HW3_APD

## Implementation Overview

This project demonstrates the usefulness of MPI (Message Passing Interface) in distributed programming. The implementation focuses on creating a topology of leader and worker processes, distributing work, and aggregating results.

### Key Components

1. **Topology Creation**: 
   - Leaders obtain their worker processes from an input file.
   - A vector of vectors represents the topology.
   - Each Cluster object maintains its own workers in a private vector.

2. **Communication Setup**:
   - Leaders send their rank to worker processes for communication.
   - Leaders exchange worker process information with neighboring leaders.

3. **Vector Processing**:
   - Leader 0 generates the initial vector and distributes it:
     - To its own worker processes.
     - To the right-side leader.
   - This process is repeated among leaders.
   - The goal is to obtain the calculated vector in process 0.

4. **Worker Process Operations**:
   - Receive the vector and process it based on the given relationship.
   - Pass the processed vector back to their leader.

5. **Result Aggregation**:
   - Intermediate processes transfer obtained vectors to the left side.
   - Leader 0 displays the final vector using the MPI 'recv()' function.

### Topology and Vector Calculation

- The approach for studying topology and vector calculation is similar.
- Key difference: The traversal direction of the graph/topology.
  - Leader 0: Can only send/receive from the right side.
  - Leaders 2 and 3: Can communicate from both sides.
  - Leader 1: Can only communicate from the left side.
- Vector route changes: 0 -> 1 -> 2 -> 3 becomes 0 -> 3 -> 2 -> 1.

## Bonus Implementation

The bonus implementation attempts to isolate process 1 from others:

- Process 1 can only communicate with its own worker processes.
- Ignoring leader 1 also excludes its children from the main computation.
- For participating processes:
  - Calculate the current process rank and number of workers.
  - Use the same relationship from the problem statement.
- Determine these by traversing the topology studied by each participating process.

## References

- [Distributed Programming Examples (MolSSI)](https://education.molssi.org/parallel-programming/04-distributed-examples/index.html)
- [Open MPI FAQ: Troubleshooting](https://www.open-mpi.org/faq/?category=troubleshooting)
- [Introduction to MPI](https://rantahar.github.io/introduction-to-mpi/aio/index.html)
- [C++ Forum: MPI Discussion](https://cplusplus.com/forum/general/282234/)
