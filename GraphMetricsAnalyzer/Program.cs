using System;
using System.Collections.Generic;
using System.Linq;
using QuickGraph;

namespace SocialNetworkGraphAnalysis
{
    class Program
    {
        static void Main(string[] args)
        {
            // Задаем параметры случайного графа
            int vertexCount = 10; // Количество вершин
            double edgeProbability = 0.3; // Вероятность появления ребра между двумя вершинами

            // Генерируем случайный граф
            var graph = GenerateRandomGraph(vertexCount, edgeProbability);

            // Выводим информацию о графе
            Console.WriteLine($"Случайный граф с {vertexCount} вершинами и вероятностью ребра {edgeProbability}:");
            Console.WriteLine($"Количество ребер: {graph.EdgeCount}");

            // Вычисляем метрики
            foreach (var vertex in graph.Vertices)
            {
                Console.WriteLine($"Вершина {vertex}:");
                Console.WriteLine($"  Степень: {GetDegree(graph, vertex)}");
                Console.WriteLine($"  Взвешенная степень: {GetWeightedDegree(graph, vertex)}");
                Console.WriteLine($"  Коэффициент кластеризации: {GetClusteringCoefficient(graph, vertex):F2}");
            }

            Console.WriteLine($"Коэффициент ассортативности: {GetAssortativityCoefficient(graph):F2}");
            Console.WriteLine($"Плотность графа: {GetGraphDensity(graph):F2}");

            // Вычисляем модулярность
            var modularity = GetModularity(graph);
            Console.WriteLine($"Модулярность: {modularity:F2}");
        }

        // Генерация случайного графа
        static UndirectedGraph<int, Edge<int>> GenerateRandomGraph(int vertexCount, double edgeProbability)
        {
            var graph = new UndirectedGraph<int, Edge<int>>();
            var random = new Random();

            // Добавляем вершины
            for (int i = 1; i <= vertexCount; i++)
            {
                graph.AddVertex(i);
            }

            // Добавляем ребра случайным образом
            for (int i = 1; i <= vertexCount; i++)
            {
                for (int j = i + 1; j <= vertexCount; j++)
                {
                    if (random.NextDouble() < edgeProbability)
                    {
                        graph.AddEdge(new Edge<int>(i, j));
                    }
                }
            }

            return graph;
        }

        // Степень узла
        static int GetDegree(UndirectedGraph<int, Edge<int>> graph, int vertex)
        {
            return graph.AdjacentEdges(vertex).Count();
        }

        // Взвешенная степень узла (в данном случае просто степень, так как веса не заданы)
        static int GetWeightedDegree(UndirectedGraph<int, Edge<int>> graph, int vertex)
        {
            return graph.AdjacentEdges(vertex).Count();
        }

        // Коэффициент кластеризации узла
        static double GetClusteringCoefficient(UndirectedGraph<int, Edge<int>> graph, int vertex)
        {
            var neighbors = GetNeighbors(graph, vertex).ToList();
            if (neighbors.Count < 2) return 0; // Если соседей меньше двух, коэффициент равен 0

            // Подсчитываем количество ребер между соседями
            int connectedPairs = 0;
            for (int i = 0; i < neighbors.Count; i++)
            {
                for (int j = i + 1; j < neighbors.Count; j++)
                {
                    if (graph.ContainsEdge(neighbors[i], neighbors[j]))
                    {
                        connectedPairs++;
                    }
                }
            }

            // Формула коэффициента кластеризации
            return (double)connectedPairs / (neighbors.Count * (neighbors.Count - 1) / 2);
        }

        // Получение соседей узла
        static IEnumerable<int> GetNeighbors(UndirectedGraph<int, Edge<int>> graph, int vertex)
        {
            var neighbors = new HashSet<int>();
            foreach (var edge in graph.AdjacentEdges(vertex))
            {
                if (edge.Source == vertex)
                    neighbors.Add(edge.Target);
                else if (edge.Target == vertex)
                    neighbors.Add(edge.Source);
            }
            return neighbors;
        }

        // Коэффициент ассортативности
        static double GetAssortativityCoefficient(UndirectedGraph<int, Edge<int>> graph)
        {
            var degrees = graph.Vertices.ToDictionary(v => v, v => GetDegree(graph, v));
            double numerator = 0, denominator1 = 0, denominator2 = 0;

            foreach (var edge in graph.Edges)
            {
                int d1 = degrees[edge.Source];
                int d2 = degrees[edge.Target];
                numerator += (d1 * d2);
                denominator1 += (d1 * d1);
                denominator2 += (d2 * d2);
            }

            double denominator = Math.Sqrt(denominator1 * denominator2);
            return numerator / denominator;
        }

        // Плотность графа
        static double GetGraphDensity(UndirectedGraph<int, Edge<int>> graph)
        {
            int vertexCount = graph.VertexCount;
            int edgeCount = graph.EdgeCount;
            return (double)(2 * edgeCount) / (vertexCount * (vertexCount - 1));
        }

        // Модулярность
        static double GetModularity(UndirectedGraph<int, Edge<int>> graph)
        {
            // Разбиваем граф на сообщества с помощью алгоритма Лувена
            var communities = LouvainCommunityDetection(graph);

            // Вычисляем модулярность
            double m = graph.EdgeCount;
            double modularity = 0;

            foreach (var community in communities)
            {
                var communityVertices = community.ToList();
                for (int i = 0; i < communityVertices.Count; i++)
                {
                    for (int j = i + 1; j < communityVertices.Count; j++)
                    {
                        int v1 = communityVertices[i];
                        int v2 = communityVertices[j];
                        int k1 = GetDegree(graph, v1);
                        int k2 = GetDegree(graph, v2);

                        double Aij = graph.ContainsEdge(v1, v2) ? 1 : 0;
                        modularity += (Aij - (k1 * k2) / (2 * m)) / (2 * m);
                    }
                }
            }

            return modularity;
        }

        // Алгоритм Лувена для разбиения графа на сообщества
        static List<HashSet<int>> LouvainCommunityDetection(UndirectedGraph<int, Edge<int>> graph)
        {
            // Инициализация: каждая вершина в своем сообществе
            var communities = graph.Vertices.Select(v => new HashSet<int> { v }).ToList();

            bool changed;
            do
            {
                changed = false;

                // Перемещаем вершины между сообществами
                foreach (var vertex in graph.Vertices)
                {
                    var currentCommunity = communities.First(c => c.Contains(vertex));
                    var bestCommunity = currentCommunity;
                    double bestGain = 0;

                    // Пробуем переместить вершину в соседние сообщества
                    foreach (var neighbor in GetNeighbors(graph, vertex))
                    {
                        var neighborCommunity = communities.First(c => c.Contains(neighbor));
                        if (neighborCommunity == currentCommunity) continue;

                        // Вычисляем прирост модулярности
                        double gain = CalculateModularityGain(graph, vertex, currentCommunity, neighborCommunity);
                        if (gain > bestGain)
                        {
                            bestGain = gain;
                            bestCommunity = neighborCommunity;
                        }
                    }

                    // Перемещаем вершину в лучшее сообщество
                    if (bestCommunity != currentCommunity)
                    {
                        currentCommunity.Remove(vertex);
                        bestCommunity.Add(vertex);
                        changed = true;
                    }
                }

                // Объединяем сообщества, если это улучшает модулярность
                communities = MergeCommunities(communities);

            } while (changed);

            return communities;
        }

        // Вычисление прироста модулярности при перемещении вершины
        static double CalculateModularityGain(UndirectedGraph<int, Edge<int>> graph, int vertex, HashSet<int> currentCommunity, HashSet<int> targetCommunity)
        {
            double m = graph.EdgeCount;
            int k_i = GetDegree(graph, vertex);
            int sigma_tot = targetCommunity.Sum(v => GetDegree(graph, v));
            int sigma_in = targetCommunity.Sum(v => graph.ContainsEdge(vertex, v) ? 1 : 0);

            return (sigma_in / (2.0 * m)) - (sigma_tot * k_i) / (2 * m * m);
        }

        // Объединение сообществ
        static List<HashSet<int>> MergeCommunities(List<HashSet<int>> communities)
        {
            var mergedCommunities = new List<HashSet<int>>();
            foreach (var community in communities)
            {
                if (community.Count > 0)
                {
                    mergedCommunities.Add(community);
                }
            }
            return mergedCommunities;
        }
    }
}