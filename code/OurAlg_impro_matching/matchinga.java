package OurAlg_impro_matching;

import Tools.*;

import java.util.*;
import java.util.List;

import static Tools.tools.distance;

public class matchinga {
    private static List<Integer> mainCenter, thisCenter, thisCL;
    private static final int NIL = 0;
    public static List<Integer>[] edgesList;
    private static final int INF = Integer.MAX_VALUE;
    public static Set<Integer> Y_prime;
    public static Set<Integer> C_prime;
    private static int[] pair, edgeSize, dist;
    private static int uN, vN;
    private static Set<Integer> Ynew;


    private static boolean bfs() {
        Queue<Integer> queue = new ArrayDeque<>();
        for (int u = 1; u <= uN; u++) {
            if (pair[u] == NIL) {
                dist[u] = 0;
                queue.add(u);
            } else {
                dist[u] = INF;
            }
        }
        dist[NIL] = INF;

        while (!queue.isEmpty()) {
            int u = queue.poll();
            if (dist[u] < dist[NIL]) {
                for (int v : edgesList[u]) {
                    if (dist[pair[v]] == INF) {
                        dist[pair[v]] = dist[u] + 1;
                        queue.add(pair[v]);
                    }
                }
            }
        }
        return dist[NIL] != INF;
    }

    private static boolean dfs(int u) {
        if (u != NIL) {
            for (int v : edgesList[u]) {
                if (dist[pair[v]] == dist[u] + 1 && dfs(pair[v])) {
                    pair[u] = v;
                    pair[v] = u;
                    return true;
                }
            }
            dist[u] = INF;
            return false;
        }
        return true;

    }

    //Algorithm
    public static boolean maxMatching() {
        pair = new int[uN + vN + 1];
        dist = new int[uN + vN + 1];

        Y_prime = new HashSet<>();
        C_prime = new HashSet<>();


        Set<Integer> Y0 = new HashSet<>();
        int matching = 0;

        while (bfs()) {
            for (int u = 1; u <= uN; u++) {
                if (pair[u] == NIL && dfs(u)) {
                    matching++;
                }
            }
        }


        if (matching < vN) {
            for (int v = 1; v <= vN; v++) {
                if (pair[v + uN] == NIL) {
                    Y0.add(v + uN);
                }
            }

            f(Y0);

            return true;
        }
        return false;
    }

    private static void f(Set<Integer> Yi) {
        Y_prime.addAll(Yi);
        Ynew = new HashSet<>();
        for (int yi : Yi) {
            thisCenter.add(thisCL.get(yi - uN - 1));
            C_prime.addAll(edgesList[yi]);
            for (int j : edgesList[yi]) {
                thisCenter.remove(mainCenter.get((j - 1)));
                if (pair[j] != 0)
                    Ynew.add(pair[j]);
            }
        }
        if (!Y_prime.containsAll(Ynew)) {
            Ynew.removeAll(Y_prime);
            f(Ynew);
        }
    }


    public static boolean input(ArrayList<Point> pointList, ArrayList<Integer> center, List<Integer> CL, ArrayList<ArrayList<Integer>> mustLinks, float R) {


        mainCenter = new ArrayList<>(center);
        thisCenter = center;
        thisCL = CL;

        uN = thisCenter.size();
        vN = thisCL.size();
        edgeSize = new int[vN + uN + 1];
        edgesList = new List[vN + uN + 1];
        Arrays.setAll(edgesList, i -> new ArrayList<>());


        for (int i = 0; i < uN; i++) {
            Point centerPoint = pointList.get(thisCenter.get(i));
            int centerMustID = centerPoint.getMustID();

            for (int j = 0; j < vN; j++) {

                Point clPoint = pointList.get(thisCL.get(j));
                int clMustID = clPoint.getMustID();


                float maxDist = 0;
                if (clMustID != -1) {
                    for (int linkedPointID : mustLinks.get(clMustID)) {
                        float tempDist = distance(centerPoint, pointList.get(linkedPointID));
                        if (centerMustID != -1) {
                            int k = 0;
                            for (int mustID : mustLinks.get(centerMustID)) {
                                float dismax = distance(pointList.get(linkedPointID), pointList.get(mustID));
                                tempDist = (tempDist >= dismax) ? tempDist : dismax;
                            }
                        }
//                        maxDist = Math.max(maxDist, tempDist);
                        maxDist = (maxDist >= tempDist) ? maxDist : tempDist;
                    }

                } else {
                    maxDist = distance(centerPoint, clPoint);
                    if (centerMustID != -1) {
                        for (int mustID : mustLinks.get(centerMustID)) {
                            maxDist = Math.max(maxDist, distance(clPoint, pointList.get(mustID)));
                        }
                    }
                }
                if (maxDist <= R) {
                    // "+1" since edgesList[0] is the special vertex
                    int from = i + 1;
                    int to = uN + j + 1;

                    edgesList[from].add(to);
                    edgesList[to].add(from);
                    edgeSize[from]++;
                    edgeSize[to]++;

                }
            }
        }
        long startTime = System.nanoTime();
        boolean a = maxMatching();
        Approx_improv_matching_kCenter.time3 += System.nanoTime() - startTime;
        return a;
    }


}
