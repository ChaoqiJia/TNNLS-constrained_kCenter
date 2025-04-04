package OurAlg_greedy_matching;

import OurAlg_LP.ApproxkCenter;
import Tools.*;

import java.util.*;

import static Tools.tools.distance;


/*
 *    HopcroftKarp algorithm
 *
 *  An algorithm that takes as input a bipartite graph
 *  and produces as output a maximum cardinality matching
 *  a set of as many edges as possible with the property
 *  that no two edges share an endpoint.
 *
 */
public class HopcroftKarp {

    private static final int nil = 0;
    private static final int infinity = Integer.MAX_VALUE;
    public static List<Integer>[] edgesList;
    public static int matching;
    private static int[] pair, edgeSize;
    private static int[] dist;
    private static int cv;
    private static int cu;
    private static List<Integer> vertexCover = new ArrayList<Integer>(), thisCenter, thisCL;
    public static HashSet<Integer> vertexCoverC, vertexCoverV;


    /*
     *  Returns true if there is an augmenting path
     */
    private static boolean BFS() {
        Queue<Integer> queue = new ArrayDeque<>();

        // First layer of vertices (set distance to 0)
        for (int v = 1; v <= cv; v++) {
            // If this is a free vertex, add it to queue
            if (pair[v] == nil) {
                dist[v] = 0;  // v is not matched
                queue.add(v);
            } else {
                dist[v] = infinity;  // infinity so that this vertex is considered next time
            }
        }

        // Initialize distance to NIL as infinity
        dist[nil] = infinity;

        // Q is going to contain vertices of the left side only
        while (!queue.isEmpty()) {
            int v = queue.poll();

            // If this node is not NIL and can provide a shorter path to NIL
            if (dist[v] < dist[nil]) {

                for (int u : edgesList[v]) {

                    // If pair of u is not considered so far
                    // (u, pair[u]) is not yet explored edge
                    if (dist[pair[u]] == infinity) {

                        // Consider the pair and add it to queue
                        dist[pair[u]] = dist[v] + 1;
                        queue.add(pair[u]);
                    }
                }
            }
        }

        // If we could come back to NIL using alternating path of distinct vertices
        // then there is an augmenting path
        return dist[nil] != infinity;
    }


    /*
     *  Returns true if there is an augmenting path beginning with free vertex v
     */
    private static boolean DFS(int v) {
        if (v != nil) {
            for (int u : edgesList[v]) {
                // Follow the distances set by BFS
                if (dist[pair[u]] == dist[v] + 1) {
                    // If dfs for pair of u also returns true
                    if (DFS(pair[u])) {
                        pair[u] = v;
                        pair[v] = u;
                        return true;
                    }
                }
            }

            // If there is no augmenting path beginning with v
            dist[v] = infinity;
            return false;
        }
        return true;
    }

    /*
     *  Returns the size of maximum matching
     */
    private static boolean HopcroftKarp() {
        // pair[v] stores pair of v in matching. If v doesn't have any pair, then pair[v] is 0
//        boolean flag = false;
        pair = new int[cv + cu + 1];

        // dist[v] stores distance of vertices
        // dist[v] is one more than dist[v'] if v is next to v' in augmenting path
        dist = new int[cv + cu + 1];
        matching = 0;

        // Keep updating the result while there is an augmenting path
        while (BFS()) {
            // Find a free vertex
            for (int v = 1; v <= cv; v++) {
                // If current vertex is free and there is an augmenting path from current vertex
                /*
                @@@@* NEED CHECK IT
                * */
                if (pair[v] == nil && DFS(v)) {
                    matching++;
                }
            }
        }

        if (matching < cu) {
            vertexCoverC = new HashSet<>();
            vertexCoverV = new HashSet<>();
            for (int u = 1; u < cu + 1; u++) {
                if (pair[u + cv] == nil) {
                    vertexCoverV.add(u + cv);
                }
            }
            HashSet<Integer> temp = new HashSet<>(DFSFindC(vertexCoverV));
            while (!vertexCoverC.equals(temp)) {
                vertexCoverC.addAll(temp);
                vertexCoverV.addAll(DFSFindY(vertexCoverC));
                temp.addAll(DFSFindC(vertexCoverV));
            }
            temp = new HashSet<>();
//            if (vertexCoverC.size() < vertexCoverV.size()) {
//                for (int c : vertexCoverC) {
//                    vertexCover.add(thisCenter.get(c - 1));
//                }
                for (int c : vertexCoverC) {
                    temp.add(thisCenter.get(c - 1));
//                    if (vertexCoverC.contains(thisCenter.get(c - 1))) {
////                        Point point = pointList.get(thisCenter.get(c - 1));
//                        thisCenter.remove(Integer.valueOf(thisCenter.get(c - 1)));
////                        int mustID = point.getMustID();
////                        if (mustID != -1) {
////                            must.remove(Integer.valueOf(mustID));
////                        }
//                    }
                }
                thisCenter.removeAll(temp);
                for (Integer v : vertexCoverV) {
                    thisCenter.add(thisCL.get(v - cv - 1));
//                    int mustID = pointList.get(thisCL.get(v - cv - 1)).getMustID();
//                    if (mustID != -1) {
//                        must.add(mustID);
//                    }

                }

//            }
            return true;
        }
        return false;

//        return matching;
    }

    public static boolean input(ArrayList<Point> pointList, ArrayList<Integer> center, List<Integer> CL, ArrayList<ArrayList<Integer>> ML, float R) {

        thisCenter = center;
        thisCL = CL;
        cv = thisCenter.size();
        cu = thisCL.size();
        vertexCover = new ArrayList<Integer>();
        edgeSize = new int[cv + cu + 1];
        edgesList = new List[cv + cu + 1];
        Arrays.setAll(edgesList, i -> new ArrayList<>());

        for (int i = 0; i < cv; i++) {
            Point centerPoint = pointList.get(thisCenter.get(i));
            int centerMustID = centerPoint.getMustID();

            for (int j = 0; j < cu; j++) {
                //input: data
                float maxdist = 0;
                Point clPoint = pointList.get(thisCL.get(j));
                int clMustID = clPoint.getMustID();

                if (clMustID != -1) {
                    for (int ml_p : ML.get(clMustID)) {
                        float tempDist = distance(centerPoint, pointList.get(ml_p));
                        if (centerMustID != -1) {
                            for (int mustID : ML.get(centerMustID)) {
                                tempDist = Math.max(tempDist, distance(pointList.get(ml_p), pointList.get(mustID)));
                            }
                        }
                        maxdist = Math.max(maxdist, tempDist);
                    }
                } else {
                    maxdist = distance(centerPoint, pointList.get(thisCL.get(j)));
                    if (centerMustID != -1) {
                        for (int mustID : ML.get(centerMustID)) {
                            maxdist = Math.max(maxdist, distance(clPoint, pointList.get(mustID)));
                        }
                    }
                }


                if (maxdist <= R) {
                    // "+1" since edgesList[0] is the special vertex
                    int from = i + 1;
                    int to = j + 1;

                    edgesList[from].add(cv + to);
                    edgesList[cv + to].add(from);
                    edgeSize[from]++;
                    edgeSize[cv + to]++;

                }


            }
        }
        long startTime = System.nanoTime();
        boolean a = HopcroftKarp();
        Approx_greedy_matching_kCenter.timeGM3 += System.nanoTime() - startTime;


        return a;
    }

//    public static float distance(Point x, Point y) {
//        double dist = 0;
//        for (int i = 0; x != null && i < x.getFeature().length; i++) {
//            dist += (x.getFeature()[i] - y.getFeature()[i]) * (x.getFeature()[i] - y.getFeature()[i]);
//        }
//        return (float) Math.sqrt(dist);
//    }

    private static HashSet<Integer> DFSFindC(HashSet<Integer> v) {
        HashSet<Integer> C_prime = new HashSet<>();
        if (!v.isEmpty()) {
            for (int y : v) {
                C_prime.addAll(edgesList[y]);
            }
        }
        return C_prime;
    }

    private static HashSet<Integer> DFSFindY(HashSet<Integer> C) {
        HashSet<Integer> Y_prime = new HashSet<>();
        if (!C.isEmpty()) {
            for (int c : C) {
                if (pair[c] != 0)
                    Y_prime.add(pair[c]);
            }
        }
        return Y_prime;
    }
}
