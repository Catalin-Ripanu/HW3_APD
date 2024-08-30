#include <mpi.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <bits/stdc++.h>
using namespace std;

// O clasa utila pentru a defini un cluster

class Cluster
{
public:
    void setClusterTopology(char *fileName);
    int getRank();
    void setRank(int rank);
    vector<int> getClusterWorkers();
    void setClusterWorkers(vector<int> clusterWorkers);

    // Proprietatile unui cluster
private:
    int rank;
    vector<int> clusterWorkers;
};

// Metode simple de prelucrare si obtinere

int Cluster::getRank()
{
    return rank;
}

void Cluster::setRank(int rank)
{
    this->rank = rank;
}

void Cluster::setClusterWorkers(vector<int> clusterWorkers)
{
    this->clusterWorkers = clusterWorkers;
}

vector<int> Cluster::getClusterWorkers()
{
    return clusterWorkers;
}

// O metoda care adauga workers in vectorul cluster-ului

void Cluster::setClusterTopology(char *fileName)
{
    int nrWorkers;
    ifstream fin(fileName);

    fin >> nrWorkers;
    clusterWorkers.resize(nrWorkers);

    for (int i = 0; i < nrWorkers; i++)
    {
        fin >> clusterWorkers[i];
        MPI_Send(&rank, 1, MPI_INT, clusterWorkers[i], 0, MPI_COMM_WORLD);
        cout << "M(" << rank << "," << clusterWorkers[i] << ")\n";
    }
}

// Metode utile pentru evitarea repetarii codului

void getClustersElem(vector<vector<int>> clusters, int iter1, int iter2)
{
    if (iter1 == clusters[iter2].size() - 1)
        cout << clusters[iter2][iter1] << " ";
    else
        cout << clusters[iter2][iter1] << ",";
}

// Metoda care afiseaza topologia din perspectiva unui proces

void printFinalFormat(int rank, vector<vector<int>> clusters)
{
    cout << rank << " -> ";
    for (int k = 0; k < clusters.size(); k++)
    {
        cout << k << ":";
        for (int l = 0; l < clusters[k].size(); l++)
        {
            getClustersElem(clusters, l, k);
        }
    }
    cout << "\n";
}

// Metoda care afiseaza topologia din perspectiva unui proces (in cazul unei partitii existente)

void printFinalPartitionFormat(int rank, vector<vector<int>> clusters)
{
    cout << rank << " -> ";
    for (int k = 0; k < clusters.size(); k++)
    {
        if (clusters[k].size() != 0)
        {
            cout << k << ":";
            for (int l = 0; l < clusters[k].size(); l++)
            {
                if (l == clusters[k].size() - 1)
                    cout << clusters[k][l] << " ";
                else
                    cout << clusters[k][l] << ",";
            }
        }
    }
    cout << "\n";
}

void getBeginEnd(int &begin, int &end, int dim, int num1, int num2)
{
    begin = num1 * (double)dim / num2;
    end = fmin((num1 + 1) * (double)dim / num2, dim);
}

// Metode care afiseaza corespunzator comunicatiile intre procese

void showFormat(int rank1, int rank2)
{
    cout << "M(" << rank1 << "," << rank2 << ")"
         << "\n";
}

void calcPartitionCluster(int &begin, int &end, int dim, int num1, int num2, vector<vector<int>> clusters)
{
    int cont = num2;
    int aux = 0;
    for (int l = 0; l < num1; l++)
    {
        if (l != 1)
        {
            cont += clusters[l].size();
        }
    }
    for (int l = 0; l < clusters.size(); l++)
    {
        aux += clusters[l].size();
    }

    getBeginEnd(begin, end, dim, cont, aux);
}

// Metode utile pentru a studia problema in diferite scenarii

class Topology
{
public:
    static void noRank0ErrorTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status);

    static void rank0ErrorTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status);

    static void noRank3ErrorTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status);

    static void rank1ErrorTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status);

    static void rank3ErrorTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status);

    static void rank2ErrorTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status);

    static void noErrorTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status);

    static void workerTopology(vector<vector<int>> &clusters, int rank, int leader, MPI_Status &status);

    static void calc0Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                              vector<int> &result, int dim, int numberProc);

    static void calcTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                             vector<int> &result, int dim, int numberProc);

    static void calc3Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                              vector<int> &result, int dim, int numberProc);

    static void calcWorkerTopology(vector<vector<int>> &clusters, int rank, MPI_Status &status,
                                   vector<int> &result, int dim, int numberProc, int leader);

    static void errorCalc0Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                                   vector<int> &result, int dim, int numberProc);

    static void errorCalc1Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                                   vector<int> &result, int dim, int numberProc);

    static void errorCalc2Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                                   vector<int> &result, int dim, int numberProc);

    static void errorCalc3Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                                   vector<int> &result, int dim, int numberProc);

    static void partitionRank0Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status);

    static void partitionRank1Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status);

    static void partitionRank2Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status);

    static void partitionRank3Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status);

    static void partitionWorkerTopology(vector<vector<int>> &clusters, int rank, int leader, MPI_Status &status);

    static void partitionCalc0Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                                       vector<int> &result, int dim, int numberProc);

    static void partitionCalc2Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                                       vector<int> &result, int dim, int numberProc);

    static void partitionCalc3Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                                       vector<int> &result, int dim, int numberProc);

    static void partitionCalcWorkersTopology(vector<vector<int>> &clusters, int rank, MPI_Status &status,
                                             vector<int> &result, int dim, int numberProc, int leader);
};

void calcSpecialPartitionCluster(int &begin, int &end, int dim, int rank, int leader, vector<vector<int>> clusters)
{
    auto findRank = find(clusters[leader].begin(), clusters[leader].end(), rank);
    int aux = 0;
    aux = aux + (int)(findRank - clusters[leader].begin());
    for (int k = 0; k < leader; k++)
    {
        if (k != 1)
        {
            aux = aux + clusters[k].size();
        }
    }
    int var = 0;
    for (int l = 0; l < clusters.size(); l++)
    {
        var = var + clusters[l].size();
    }
    getBeginEnd(begin, end, dim, aux, var);
}

// Metode care prelucreaza un cluster dat

void prepCluster(Cluster *cluster, int rank, vector<int> &result,
                 int dim, MPI_Status &status, int numberProc)
{
    for (int i = 0; i < cluster->getClusterWorkers().size(); i++)
    {
        MPI_Send(&result[0], dim, MPI_INT, cluster->getClusterWorkers()[i], rank, MPI_COMM_WORLD);
        showFormat(cluster->getRank(), cluster->getClusterWorkers()[i]);
    }
    for (int i = 0; i < cluster->getClusterWorkers().size(); i++)
    {
        int begin, end;
        getBeginEnd(begin, end, dim, cluster->getClusterWorkers()[i] - 4, numberProc - 4);
        MPI_Recv(&result[begin], end - begin, MPI_INT, cluster->getClusterWorkers()[i],
                 cluster->getClusterWorkers()[i], MPI_COMM_WORLD, &status);
    }
}

void calcCluster0(int iterator, vector<int> &helper, int val,
                  Cluster *cluster, vector<vector<int>> &clusters, MPI_Status &status, int rank1, int rank2)
{
    if (iterator != cluster->getRank())
    {
        MPI_Recv(&helper[0], val, MPI_INT, rank1, iterator, MPI_COMM_WORLD, &status);
        helper.resize(helper.size() - count(helper.begin(), helper.end(), 0));
        clusters[iterator] = helper;
        MPI_Send(&helper[0], helper.size(), MPI_INT, rank2, iterator, MPI_COMM_WORLD);
        showFormat(cluster->getRank(), rank2);
        for (int j = 0; j < cluster->getClusterWorkers().size(); j++)
        {
            MPI_Send(&helper[0], helper.size(), MPI_INT, cluster->getClusterWorkers()[j], iterator, MPI_COMM_WORLD);
            showFormat(cluster->getRank(), cluster->getClusterWorkers()[j]);
        }
    }
}

void workCluster0(Cluster *cluster, int paramRank1, int paramRank2, int paramRank3, int paramRank4)
{
    MPI_Send(&cluster->getClusterWorkers()[0],
             cluster->getClusterWorkers().size(), MPI_INT,
             paramRank1, paramRank2, MPI_COMM_WORLD);
    showFormat(cluster->getRank(), paramRank3);
    for (int i = 0; i < cluster->getClusterWorkers().size(); i++)
    {
        MPI_Send(&cluster->getClusterWorkers()[0],
                 cluster->getClusterWorkers().size(), MPI_INT,
                 cluster->getClusterWorkers()[i], paramRank4, MPI_COMM_WORLD);
        showFormat(cluster->getRank(), cluster->getClusterWorkers()[i]);
    }
}

// Calculul vectorului (cand exista o partitie) de catre procesele workers

void Topology::partitionCalcWorkersTopology(vector<vector<int>> &clusters, int rank, MPI_Status &status,
                                            vector<int> &result, int dim, int numberProc, int leader)
{
    if (leader != 1 && rank != 1)
    {
        MPI_Recv(&result[0], dim, MPI_INT, leader, leader, MPI_COMM_WORLD, &status);
        int begin, end;
        calcSpecialPartitionCluster(begin, end, dim, rank, leader, clusters);
        for (int i = begin; i < end; i++)
        {
            result[i] = result[i] * 5;
        }
        MPI_Send(&result[begin], end - begin, MPI_INT, leader, rank, MPI_COMM_WORLD);
        showFormat(rank, leader);
    }
}

// Calculul vectorului (cand exista o partitie) de catre procesul 3

void Topology::partitionCalc3Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                                      vector<int> &result, int dim, int numberProc)
{
    if (cluster->getRank() == 3)
    {
        MPI_Recv(&result[0], dim, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

        MPI_Send(&result[0], dim, MPI_INT, cluster->getRank() - 1, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), cluster->getRank() - 1);

        for (int i = 0; i < cluster->getClusterWorkers().size(); i++)
        {
            MPI_Send(&result[0], dim, MPI_INT, cluster->getClusterWorkers()[i], cluster->getRank(), MPI_COMM_WORLD);
            showFormat(cluster->getRank(), cluster->getClusterWorkers()[i]);
        }
        for (int i = 0; i < cluster->getClusterWorkers().size(); i++)
        {
            int begin, end;
            calcPartitionCluster(begin, end, dim, cluster->getRank(), i, clusters);
            MPI_Recv(&result[begin], end - begin, MPI_INT, cluster->getClusterWorkers()[i],
                     cluster->getClusterWorkers()[i], MPI_COMM_WORLD, &status);
        }
        MPI_Send(&result[0], dim, MPI_INT, 0, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), 0);
        MPI_Recv(&result[0], dim, MPI_INT, cluster->getRank() - 1, cluster->getRank() - 1, MPI_COMM_WORLD, &status);
        MPI_Send(&result[0], dim, MPI_INT, 0, cluster->getRank() - 1, MPI_COMM_WORLD);
        showFormat(cluster->getRank(), 0);
    }
}

// Calculul vectorului (cand exista o partitie) de catre procesul 2

void Topology::partitionCalc2Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                                      vector<int> &result, int dim, int numberProc)
{
    if (cluster->getRank() == 2)
    {
        MPI_Recv(&result[0], dim, MPI_INT, cluster->getRank() + 1, cluster->getRank() + 1, MPI_COMM_WORLD, &status);

        for (int i = 0; i < cluster->getClusterWorkers().size(); i++)
        {
            MPI_Send(&result[0], dim, MPI_INT, cluster->getClusterWorkers()[i], cluster->getRank(), MPI_COMM_WORLD);
            showFormat(cluster->getRank(), cluster->getClusterWorkers()[i]);
        }
        for (int i = 0; i < cluster->getClusterWorkers().size(); i++)
        {
            int begin, end;
            calcPartitionCluster(begin, end, dim, cluster->getRank(), i, clusters);
            MPI_Recv(&result[begin], end - begin, MPI_INT, cluster->getClusterWorkers()[i],
                     cluster->getClusterWorkers()[i], MPI_COMM_WORLD, &status);
        }
        MPI_Send(&result[0], dim, MPI_INT, cluster->getRank() + 1, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), cluster->getRank() + 1);
    }
}

// Calculul vectorului (cand exista o partitie) de catre procesul 0

void Topology::partitionCalc0Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                                      vector<int> &result, int dim, int numberProc)
{
    if (cluster->getRank() == 0)
    {
        for (int i = 0; i < result.size(); i++)
        {
            result[i] = result.size() - i - 1;
        }
        MPI_Send(&result[0], dim, MPI_INT, 3, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), 3);
        for (int i = 0; i < cluster->getClusterWorkers().size(); i++)
        {
            MPI_Send(&result[0], dim, MPI_INT, cluster->getClusterWorkers()[i], cluster->getRank(), MPI_COMM_WORLD);
            showFormat(cluster->getRank(), cluster->getClusterWorkers()[i]);
        }
        for (int i = 0; i < cluster->getClusterWorkers().size(); i++)
        {

            int begin, end;
            calcPartitionCluster(begin, end, dim, i, cluster->getRank(), clusters);
            MPI_Recv(&result[begin], end - begin, MPI_INT, cluster->getClusterWorkers()[i],
                     cluster->getClusterWorkers()[i], MPI_COMM_WORLD, &status);
        }
        for (int i = 3; i > 1; i--)
        {
            vector<int> helper(dim, 0);
            MPI_Recv(&helper[0], dim, MPI_INT, 3, i, MPI_COMM_WORLD, &status);
            for (int k = 0; k < clusters[i].size(); k++)
            {
                int begin, end;
                calcPartitionCluster(begin, end, dim, i, k, clusters);
                for (int l = begin; l < end; l++)
                {
                    result[l] = helper[l];
                }
            }
        }

        cout << "Rezultat: ";
        for (int l = 0; l < result.size(); l++)
        {
            cout << result[l] << " ";
        }
        cout << "\n";
    }
}

// Studiul topologiei (cand exista o partitie) de catre procesele workers

void Topology::partitionWorkerTopology(vector<vector<int>> &clusters, int rank, int leader, MPI_Status &status)
{
    int val = 10;
    if (leader == 1)
    {
        vector<int> helper(val, 0);
        MPI_Recv(&helper[0], val, MPI_INT, leader, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        helper.resize(helper.size() - count(helper.begin(), helper.end(), 0));
        clusters[status.MPI_TAG] = helper;
    }
    else
    {
        for (int i = 0; i < 3; i++)
        {
            vector<int> helper(val, 0);
            MPI_Recv(&helper[0], val, MPI_INT, leader, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            helper.resize(helper.size() - count(helper.begin(), helper.end(), 0));
            clusters[status.MPI_TAG] = helper;
        }
    }

    printFinalPartitionFormat(rank, clusters);
}

// Studiul topologiei (cand exista o partitie) de catre procesul 2

void Topology::partitionRank2Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status)
{
    if (cluster->getRank() == 2)
    {
        int val = 10;
        workCluster0(cluster, cluster->getRank() + 1, cluster->getRank(), cluster->getRank() + 1, cluster->getRank());

        for (int i = 0; i < 2; i++)
        {
            vector<int> helper(val, 0);
            MPI_Recv(&helper[0], val, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            helper.resize(helper.size() - count(helper.begin(), helper.end(), 0));
            clusters[status.MPI_TAG] = helper;
            for (int j = 0; j < cluster->getClusterWorkers().size(); j++)
            {
                MPI_Send(&helper[0], helper.size(), MPI_INT, cluster->getClusterWorkers()[j], status.MPI_TAG, MPI_COMM_WORLD);
                showFormat(cluster->getRank(), cluster->getClusterWorkers()[j]);
            }
        }

        printFinalPartitionFormat(cluster->getRank(), clusters);
    }
}

// Studiul topologiei (cand exista o partitie) de catre procesul 3

void Topology::partitionRank3Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status)
{
    if (cluster->getRank() == 3)
    {
        int val = 10;
        MPI_Send(&cluster->getClusterWorkers()[0],
                 cluster->getClusterWorkers().size(), MPI_INT,
                 cluster->getRank() - 1, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), cluster->getRank() - 1);
        workCluster0(cluster, 0, cluster->getRank(), 0, cluster->getRank());

        for (int i = 0; i < 2; i++)
        {
            vector<int> helper(val, 0);
            MPI_Recv(&helper[0], val, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            helper.resize(helper.size() - count(helper.begin(), helper.end(), 0));
            clusters[status.MPI_TAG] = helper;
            if (status.MPI_SOURCE == 0)
            {
                MPI_Send(&helper[0], helper.size(), MPI_INT, cluster->getRank() - 1, status.MPI_TAG, MPI_COMM_WORLD);
                showFormat(cluster->getRank(), cluster->getRank() - 1);
            }
            else
            {
                MPI_Send(&helper[0], helper.size(), MPI_INT, 0, status.MPI_TAG, MPI_COMM_WORLD);
                showFormat(cluster->getRank(), 0);
            }
            for (int j = 0; j < cluster->getClusterWorkers().size(); j++)
            {
                MPI_Send(&helper[0], helper.size(), MPI_INT, cluster->getClusterWorkers()[j], status.MPI_TAG, MPI_COMM_WORLD);
                showFormat(cluster->getRank(), cluster->getClusterWorkers()[j]);
            }
        }

        printFinalPartitionFormat(cluster->getRank(), clusters);
    }
}

// Studiul topologiei (cand exista o partitie) de catre procesul 1

void Topology::partitionRank1Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status)
{
    if (cluster->getRank() == 1)
    {
        for (int i = 0; i < cluster->getClusterWorkers().size(); i++)
        {
            MPI_Send(&cluster->getClusterWorkers()[0],
                     cluster->getClusterWorkers().size(), MPI_INT,
                     cluster->getClusterWorkers()[i], cluster->getRank(), MPI_COMM_WORLD);
            showFormat(cluster->getRank(), cluster->getClusterWorkers()[i]);
        }

        printFinalPartitionFormat(cluster->getRank(), clusters);
    }
}

// Studiul topologiei (cand exista o partitie) de catre procesul 0

void Topology::partitionRank0Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status)
{
    if (cluster->getRank() == 0)
    {
        int val = 10;
        workCluster0(cluster, 3, 0, 3, 0);

        for (int i = 0; i < 2; i++)
        {
            vector<int> helper(val, 0);
            MPI_Recv(&helper[0], val, MPI_INT, 3, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            helper.resize(helper.size() - count(helper.begin(), helper.end(), 0));
            clusters[status.MPI_TAG] = helper;
            for (int j = 0; j < cluster->getClusterWorkers().size(); j++)
            {
                MPI_Send(&helper[0], helper.size(), MPI_INT, cluster->getClusterWorkers()[j], status.MPI_TAG, MPI_COMM_WORLD);
                showFormat(cluster->getRank(), cluster->getClusterWorkers()[j]);
            }
        }
        printFinalPartitionFormat(cluster->getRank(), clusters);
    }
}

// Calculul vectorului (cu erori) de catre procesul 3

void Topology::errorCalc3Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                                  vector<int> &result, int dim, int numberProc)
{
    if (cluster->getRank() == 3)
    {
        MPI_Recv(&result[0], dim, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Send(&result[0], dim, MPI_INT, cluster->getRank() - 1, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), cluster->getRank() - 1);
        prepCluster(cluster, cluster->getRank(), result, dim, status, numberProc);
        MPI_Send(&result[0], dim, MPI_INT, 0, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), 0);
        for (int i = 2; i > 0; i--)
        {
            MPI_Recv(&result[0], dim, MPI_INT, cluster->getRank() - 1, i, MPI_COMM_WORLD, &status);
            MPI_Send(&result[0], dim, MPI_INT, 0, i, MPI_COMM_WORLD);
            showFormat(cluster->getRank(), 0);
        }
    }
}

// Calculul vectorului (cu erori) de catre procesul 2

void Topology::errorCalc2Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                                  vector<int> &result, int dim, int numberProc)
{
    if (cluster->getRank() == 2)
    {
        MPI_Recv(&result[0], dim, MPI_INT, cluster->getRank() + 1, cluster->getRank() + 1, MPI_COMM_WORLD, &status);
        MPI_Send(&result[0], dim, MPI_INT, cluster->getRank() - 1, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), cluster->getRank() - 1);
        prepCluster(cluster, cluster->getRank(), result, dim, status, numberProc);
        MPI_Send(&result[0], dim, MPI_INT, cluster->getRank() + 1, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), cluster->getRank() + 1);
        MPI_Recv(&result[0], dim, MPI_INT, cluster->getRank() - 1, cluster->getRank() - 1, MPI_COMM_WORLD, &status);
        MPI_Send(&result[0], dim, MPI_INT, cluster->getRank() + 1, cluster->getRank() - 1, MPI_COMM_WORLD);
        showFormat(cluster->getRank(), cluster->getRank() + 1);
    }
}

// Calculul vectorului (cu erori) de catre procesul 1

void Topology::errorCalc1Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                                  vector<int> &result, int dim, int numberProc)
{
    if (cluster->getRank() == 1)
    {
        MPI_Recv(&result[0], dim, MPI_INT, cluster->getRank() + 1, cluster->getRank() + 1, MPI_COMM_WORLD, &status);
        prepCluster(cluster, cluster->getRank(), result, dim, status, numberProc);
        MPI_Send(&result[0], dim, MPI_INT, cluster->getRank() + 1, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), cluster->getRank() + 1);
    }
}

// Calculul vectorului (cu erori) de catre procesul 0

void Topology::errorCalc0Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                                  vector<int> &result, int dim, int numberProc)
{
    if (cluster->getRank() == 0)
    {
        for (int i = 0; i < result.size(); i++)
        {
            result[i] = result.size() - i - 1;
        }
        MPI_Send(&result[0], dim, MPI_INT, 3, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), 3);
        prepCluster(cluster, cluster->getRank(), result, dim, status, numberProc);
        for (int i = 3; i > 0; i--)
        {
            vector<int> helper(dim, 0);
            MPI_Recv(&helper[0], dim, MPI_INT, 3, i, MPI_COMM_WORLD, &status);
            for (int k = 0; k < clusters[i].size(); k++)
            {
                int begin, end;
                getBeginEnd(begin, end, dim, clusters[i][k] - 4, numberProc - 4);
                for (int l = begin; l < end; l++)
                {
                    result[l] = helper[l];
                }
            }
        }

        cout << "Rezultat: ";
        for (int l = 0; l < result.size(); l++)
        {
            cout << result[l] << " ";
        }
        cout << "\n";
    }
}

// Studiul topologiei (cu erori) de catre procesul 2

void Topology::rank2ErrorTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status)
{
    if (cluster->getRank() == 2)
    {
        int val = 10;
        MPI_Send(&cluster->getClusterWorkers()[0],
                 cluster->getClusterWorkers().size(), MPI_INT,
                 cluster->getRank() - 1, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), cluster->getRank() - 1);
        workCluster0(cluster, cluster->getRank() + 1,
                     cluster->getRank(), cluster->getRank() + 1, cluster->getRank());
        for (int i = 0; i < 3; i++)
        {
            vector<int> helper(val, 0);
            MPI_Recv(&helper[0], val, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            helper.resize(helper.size() - count(helper.begin(), helper.end(), 0));
            clusters[status.MPI_TAG] = helper;
            if (status.MPI_SOURCE == cluster->getRank() + 1)
            {
                MPI_Send(&helper[0], helper.size(), MPI_INT, cluster->getRank() - 1, status.MPI_TAG, MPI_COMM_WORLD);
                showFormat(cluster->getRank(), cluster->getRank() - 1);
            }
            else
            {
                MPI_Send(&helper[0], helper.size(), MPI_INT, cluster->getRank() + 1, status.MPI_TAG, MPI_COMM_WORLD);
                showFormat(cluster->getRank(), cluster->getRank() + 1);
            }
            for (int j = 0; j < cluster->getClusterWorkers().size(); j++)
            {
                MPI_Send(&helper[0], helper.size(), MPI_INT, cluster->getClusterWorkers()[j], status.MPI_TAG, MPI_COMM_WORLD);
                showFormat(cluster->getRank(), cluster->getClusterWorkers()[j]);
            }
        }
        printFinalFormat(cluster->getRank(), clusters);
    }
}

// Studiul topologiei (cu erori) de catre procesul 3

void Topology::rank3ErrorTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status)
{
    if (cluster->getRank() == 3)
    {
        int val = 10;
        MPI_Send(&cluster->getClusterWorkers()[0],
                 cluster->getClusterWorkers().size(), MPI_INT,
                 cluster->getRank() - 1, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), cluster->getRank() - 1);
        workCluster0(cluster, 0, cluster->getRank(), 0, cluster->getRank());
        for (int i = 0; i < 3; i++)
        {
            vector<int> helper(val, 0);
            MPI_Recv(&helper[0], val, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            helper.resize(helper.size() - count(helper.begin(), helper.end(), 0));
            clusters[status.MPI_TAG] = helper;
            if (status.MPI_SOURCE == 0)
            {
                MPI_Send(&helper[0], helper.size(), MPI_INT, cluster->getRank() - 1, status.MPI_TAG, MPI_COMM_WORLD);
                showFormat(cluster->getRank(), cluster->getRank() - 1);
            }
            else
            {
                MPI_Send(&helper[0], helper.size(), MPI_INT, 0, status.MPI_TAG, MPI_COMM_WORLD);
                showFormat(cluster->getRank(), 0);
            }
            for (int j = 0; j < cluster->getClusterWorkers().size(); j++)
            {
                MPI_Send(&helper[0], helper.size(), MPI_INT, cluster->getClusterWorkers()[j], status.MPI_TAG, MPI_COMM_WORLD);
                showFormat(cluster->getRank(), cluster->getClusterWorkers()[j]);
            }
        }

        printFinalFormat(cluster->getRank(), clusters);
    }
}

// Studiul topologiei (cu erori) de catre procesul 1

void Topology::rank1ErrorTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status)
{
    if (cluster->getRank() == 1)
    {
        int val = 10;
        workCluster0(cluster, cluster->getRank() + 1, cluster->getRank(), cluster->getRank() + 1, cluster->getRank());
        for (int i = 0; i < 3; i++)
        {
            vector<int> helper(val, 0);
            MPI_Recv(&helper[0], val, MPI_INT, cluster->getRank() + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            helper.resize(helper.size() - count(helper.begin(), helper.end(), 0));
            clusters[status.MPI_TAG] = helper;
            for (int j = 0; j < cluster->getClusterWorkers().size(); j++)
            {
                MPI_Send(&helper[0], helper.size(), MPI_INT, cluster->getClusterWorkers()[j], status.MPI_TAG, MPI_COMM_WORLD);
                showFormat(cluster->getRank(), cluster->getClusterWorkers()[j]);
            }
        }

        printFinalFormat(cluster->getRank(), clusters);
    }
}

// Studiul topologiei (cu erori) de catre procesul 0

void Topology::rank0ErrorTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status)
{
    if (cluster->getRank() == 0)
    {
        int val = 10;
        workCluster0(cluster, 3, 0, 3, 0);
        for (int i = 0; i < 3; i++)
        {
            vector<int> helper(val, 0);
            MPI_Recv(&helper[0], val, MPI_INT, 3, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            helper.resize(helper.size() - count(helper.begin(), helper.end(), 0));
            clusters[status.MPI_TAG] = helper;
            for (int j = 0; j < cluster->getClusterWorkers().size(); j++)
            {
                MPI_Send(&helper[0], helper.size(), MPI_INT, cluster->getClusterWorkers()[j], status.MPI_TAG, MPI_COMM_WORLD);
                showFormat(cluster->getRank(), cluster->getClusterWorkers()[j]);
            }
        }

        printFinalFormat(cluster->getRank(), clusters);
    }
}

// Calculul vectorului (fara erori / cu erori) de catre procesele workers

void Topology::calcWorkerTopology(vector<vector<int>> &clusters, int rank, MPI_Status &status,
                                  vector<int> &result, int dim, int numberProc, int leader)
{
    MPI_Recv(&result[0], dim, MPI_INT, leader, leader, MPI_COMM_WORLD, &status);
    int begin, end;
    getBeginEnd(begin, end, dim, rank - 4, numberProc - 4);
    for (int i = begin; i < end; i++)
    {
        result[i] = result[i] * 5;
    }
    MPI_Send(&result[begin], end - begin, MPI_INT, leader, rank, MPI_COMM_WORLD);
    showFormat(rank, leader);
}

// Calculul vectorului (fara erori) de catre procesul 3

void Topology::calc3Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                             vector<int> &result, int dim, int numberProc)
{
    if (cluster->getRank() == 3)
    {
        MPI_Recv(&result[0], dim, MPI_INT, cluster->getRank() - 1, cluster->getRank() - 1, MPI_COMM_WORLD, &status);
        prepCluster(cluster, cluster->getRank(), result, dim, status, numberProc);
        MPI_Send(&result[0], dim, MPI_INT, 0, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), 0);
        for (int i = 2; i > 0; i--)
        {
            MPI_Recv(&result[0], dim, MPI_INT, cluster->getRank() - 1, i, MPI_COMM_WORLD, &status);
            MPI_Send(&result[0], dim, MPI_INT, 0, i, MPI_COMM_WORLD);
            showFormat(cluster->getRank(), 0);
        }
    }
}

// Calculul vectorului (fara erori) de catre procesele 1 si 2

void Topology::calcTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                            vector<int> &result, int dim, int numberProc)
{
    if (cluster->getRank() < 3 && cluster->getRank() != 0)
    {
        MPI_Recv(&result[0], dim, MPI_INT, cluster->getRank() - 1, cluster->getRank() - 1, MPI_COMM_WORLD, &status);
        MPI_Send(&result[0], dim, MPI_INT, cluster->getRank() + 1, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), cluster->getRank() + 1);
        prepCluster(cluster, cluster->getRank(), result, dim, status, numberProc);
        MPI_Send(&result[0], dim, MPI_INT, cluster->getRank() + 1, cluster->getRank(), MPI_COMM_WORLD);
        showFormat(cluster->getRank(), cluster->getRank() + 1);
        if (cluster->getRank() == 2)
        {
            MPI_Recv(&result[0], dim, MPI_INT, cluster->getRank() - 1, cluster->getRank() - 1, MPI_COMM_WORLD, &status);
            MPI_Send(&result[0], dim, MPI_INT, cluster->getRank() + 1, cluster->getRank() - 1, MPI_COMM_WORLD);
            showFormat(cluster->getRank(), cluster->getRank() + 1);
        }
    }
}

// Calculul vectorului (fara erori) de catre procesul 0

void Topology::calc0Topology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status,
                             vector<int> &result, int dim, int numberProc)
{
    if (cluster->getRank() == 0)
    {
        for (int i = 0; i < result.size(); i++)
        {
            result[i] = result.size() - i - 1;
        }
        MPI_Send(&result[0], dim, MPI_INT, cluster->getRank() + 1, 0, MPI_COMM_WORLD);
        showFormat(cluster->getRank(), cluster->getRank() + 1);
        prepCluster(cluster, 0, result, dim, status, numberProc);
        for (int i = 3; i > 0; i--)
        {
            vector<int> helper(dim, 0);
            MPI_Recv(&helper[0], dim, MPI_INT, 3, i, MPI_COMM_WORLD, &status);
            for (int k = 0; k < clusters[i].size(); k++)
            {
                int begin, end;
                getBeginEnd(begin, end, dim, clusters[i][k] - 4, numberProc - 4);
                for (int l = begin; l < end; l++)
                {
                    result[l] = helper[l];
                }
            }
        }

        cout << "Rezultat: ";
        for (int l = 0; l < result.size(); l++)
        {
            cout << result[l] << " ";
        }
        cout << "\n";
    }
}

// Studiul topologiei (fara erori) de catre procesele workers

void Topology::workerTopology(vector<vector<int>> &clusters, int rank, int leader, MPI_Status &status)
{
    int val = 10;
    for (int i = 0; i < 4; i++)
    {
        vector<int> helper(val, 0);
        MPI_Recv(&helper[0], val, MPI_INT, leader, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        helper.resize(helper.size() - count(helper.begin(), helper.end(), 0));
        clusters[status.MPI_TAG] = helper;
    }

    printFinalFormat(rank, clusters);
}

// Studiul topologiei (fara erori) de catre procesele 1 si 2

void Topology::noErrorTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status)
{
    if (cluster->getRank() < 3 && cluster->getRank() != 0)
    {
        int val = 10;
        workCluster0(cluster, cluster->getRank() + 1, cluster->getRank(), cluster->getRank() + 1, cluster->getRank());
        for (int i = 0; i < 4; i++)
        {
            vector<int> helper(val, 0);
            calcCluster0(i, helper, val, cluster, clusters, status, cluster->getRank() - 1, cluster->getRank() + 1);
        }

        printFinalFormat(cluster->getRank(), clusters);
    }
}

// Studiul topologiei (fara erori) de catre procesul 3

void Topology::noRank3ErrorTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status)
{

    if (cluster->getRank() == 3)
    {
        int val = 10;
        workCluster0(cluster, 0, cluster->getRank(), 0, cluster->getRank());
        for (int i = 0; i < 4; i++)
        {
            vector<int> helper(val, 0);
            calcCluster0(i, helper, val, cluster, clusters, status, cluster->getRank() - 1, 0);
        }

        printFinalFormat(cluster->getRank(), clusters);
    }
}

// Studiul topologiei (fara erori) de catre procesul 0

void Topology::noRank0ErrorTopology(vector<vector<int>> &clusters, Cluster *cluster, MPI_Status &status)
{
    int val = 10;
    if (cluster->getRank() == 0)
    {
        workCluster0(cluster, cluster->getRank() + 1, 0, cluster->getRank() + 1, 0);
        for (int i = 0; i < 4; i++)
        {
            vector<int> helper(val, 0);
            calcCluster0(i, helper, val, cluster, clusters, status, 3, cluster->getRank() + 1);
        }

        printFinalFormat(cluster->getRank(), clusters);
    }
}

int main(int argc, char **argv)
{
    int numberProc, rank;
    int nrClusters = 4;
    int option;
    int dim;
    if (argv[2] == NULL)
        exit(1);
    else
    {
        dim = atoi(argv[1]);
        option = atoi(argv[2]);
    }
    vector<int> resultVector(dim, 0);
    int leader = -1;
    vector<vector<int>> vecTopology(nrClusters);
    vector<Cluster *> vecClusters(nrClusters, NULL);
    char fileName[15];
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numberProc);

    // Fiecare cluster isi obtine workers

    if (rank < 4)
    {
        Cluster *cluster = new Cluster();
        cluster->setRank(rank);
        sprintf(fileName, "./cluster%d.txt", rank);
        cluster->setClusterTopology(fileName);
        vecClusters[rank] = cluster;
        vecTopology[rank] = cluster->getClusterWorkers();
    }
    else
    {
        MPI_Recv(&leader, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
    }

    // Se asteapta procesele

    MPI_Barrier(MPI_COMM_WORLD);

    // Cazul in care nu exista erori

    if (option == 0)
    {
        // Cazul preluat de lideri (se studiaza topologia fara erori)

        if (rank < 4)
        {
            Topology::noRank0ErrorTopology(vecTopology, vecClusters[rank], status);
            Topology::noRank3ErrorTopology(vecTopology, vecClusters[rank], status);
            Topology::noErrorTopology(vecTopology, vecClusters[rank], status);
        }

        // Cazul preluat de workers

        else
        {
            Topology::workerTopology(vecTopology, rank, leader, status);
        }

        // Se asteapta procesele

        MPI_Barrier(MPI_COMM_WORLD);

        // Cazul preluat de lideri (se calculeaza acel vector, cazul fara erori)

        if (rank < 4)
        {
            Topology::calc0Topology(vecTopology, vecClusters[rank], status,
                                    resultVector, dim, numberProc);
            Topology::calcTopology(vecTopology, vecClusters[rank], status,
                                   resultVector, dim, numberProc);
            Topology::calc3Topology(vecTopology, vecClusters[rank], status,
                                    resultVector, dim, numberProc);
        }

        // Cazul preluat de workers

        else
        {
            Topology::calcWorkerTopology(vecTopology, rank, status, resultVector, dim, numberProc, leader);
        }

        // Se asteapta procesele

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Cazul in care exista erori

    else if (option == 1)
    {
        // Cazul preluat de lideri (se studiaza topologia cu erori)

        if (rank < 4)
        {
            Topology::rank0ErrorTopology(vecTopology, vecClusters[rank], status);
            Topology::rank1ErrorTopology(vecTopology, vecClusters[rank], status);
            Topology::rank3ErrorTopology(vecTopology, vecClusters[rank], status);
            Topology::rank2ErrorTopology(vecTopology, vecClusters[rank], status);
        }

        // Cazul preluat de workers

        else
        {
            Topology::workerTopology(vecTopology, rank, leader, status);
        }

        // Se asteapta procesele

        MPI_Barrier(MPI_COMM_WORLD);

        // Cazul preluat de lideri (se calculeaza acel vector, cazul cu erori in topologie)

        if (rank < 4)
        {
            Topology::errorCalc0Topology(vecTopology, vecClusters[rank], status,
                                         resultVector, dim, numberProc);
            Topology::errorCalc1Topology(vecTopology, vecClusters[rank], status,
                                         resultVector, dim, numberProc);
            Topology::errorCalc2Topology(vecTopology, vecClusters[rank], status,
                                         resultVector, dim, numberProc);
            Topology::errorCalc3Topology(vecTopology, vecClusters[rank], status,
                                         resultVector, dim, numberProc);
        }

        // Cazul preluat de workers

        else
        {
            Topology::calcWorkerTopology(vecTopology, rank, status, resultVector, dim, numberProc, leader);
        }

        // Se asteapta procesele

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Cazul in care exista o partitie

    else if (option == 2)
    {

        // Cazul preluat de lideri (se studiaza topologia existand o partitie)

        if (rank < 4)
        {
            Topology::partitionRank0Topology(vecTopology, vecClusters[rank], status);
            Topology::partitionRank1Topology(vecTopology, vecClusters[rank], status);
            Topology::partitionRank2Topology(vecTopology, vecClusters[rank], status);
            Topology::partitionRank3Topology(vecTopology, vecClusters[rank], status);
        }

        // Cazul preluat de workers

        else
        {
            Topology::partitionWorkerTopology(vecTopology, rank, leader, status);
        }

        // Se asteapta procesele

        MPI_Barrier(MPI_COMM_WORLD);

        // Cazul preluat de lideri (se calculeaza acel vector existand o partitie)

        if (rank < 4)
        {
            Topology::partitionCalc0Topology(vecTopology, vecClusters[rank], status,
                                             resultVector, dim, numberProc);
            Topology::partitionCalc2Topology(vecTopology, vecClusters[rank], status,
                                             resultVector, dim, numberProc);
            Topology::partitionCalc3Topology(vecTopology, vecClusters[rank], status,
                                             resultVector, dim, numberProc);
        }

        // Cazul preluat de workers

        else
        {
            Topology::partitionCalcWorkersTopology(vecTopology, rank, status, resultVector, dim, numberProc, leader);
        }

        // Se asteapta procesele

        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}