
package OurAlg_impro_matching;

import Tools.*;

import java.util.*;

import static Tools.tools.distance;


public class Approx_improv_matching_kCenter{

    public static int k, N;
    public static ArrayList<Point> pointList;
    public static ArrayList<ArrayList<Integer>> cannotLinkSet, mustLinkSet;
    static float[] RMark;
    public static long time3 = 0;
    static float R, errorBar, Rmax, Rmin;
    static ArrayList<Integer> must = new ArrayList<>(), centers;

    public Approx_improv_matching_kCenter(ArrayList<Point> pointList, ArrayList<ArrayList<Integer>> CannotList, ArrayList<ArrayList<Integer>> MustList, int k, float[] RMark) {
        this.pointList = pointList;
        cannotLinkSet = CannotList;
        mustLinkSet = MustList;
        this.k = k;
        this.RMark = RMark;
    }

    // Method to set clusterID for centers
    private void setClusterIDs(List<Integer> centers) {
        for (int i = 0; i < centers.size(); i++) {
            pointList.get(centers.get(i)).setClusterID(i);
        }
    }

    public float[] vertexCover(Random r) {
        System.out.println("imrds_LP:___________________________________________");
        N = pointList.size();
        float[] accuracy;
        time3 = 0;
        Rmax = RMark[0];
        Rmin = RMark[1];
        errorBar = RMark[2];
        centers = new ArrayList<>();
        must.clear();

        for (int i = 0; i < N; i++) {
            pointList.get(i).setClusterID(-1);
        }

        String times = algorithm(r);
        accuracy = Metrics.metrics(N, pointList, centers, k, times + ", " + time3);
        return accuracy;
    }

    // Method to find the candidate center set
    private List<List<Integer>> findCandidateCenterSet() {
        List<List<Integer>> candidateCenterSet = new ArrayList<>();
        int maximumSize = 0;
        for (ArrayList<Integer> integerArrayList : cannotLinkSet) {
            maximumSize = Math.max(maximumSize, integerArrayList.size());
        }
        for (ArrayList<Integer> integers : cannotLinkSet) {
            if (integers.size() == maximumSize) {
                candidateCenterSet.add(integers);
            }
        }
        return candidateCenterSet;
    }

    // Method to initialize center
    private List<Integer> initializeCenter(List<List<Integer>> candidateCenterSet, Random r) {
        List<Integer> initCenter = new ArrayList<>();
        if (candidateCenterSet.size() == 0) {
            initCenter.add(r.nextInt(N));
        } else {
            initCenter.addAll(candidateCenterSet.get(r.nextInt(candidateCenterSet.size())));
        }
        return initCenter;
    }

    // Method to perform binary search to ensure R
    private void performBinarySearch(List<Integer> initCenter) {

        boolean flag = false;
        while (Rmax - Rmin > errorBar || !flag) {
            R = (Rmin + Rmax) / 2.0f;
            if (R == Rmin || R == Rmax) {
                Rmin = R = Rmax;
//                break;
            }
            flag = true;
            centers.clear();
            centers.addAll(initCenter);

            // Exchange the center
            if (centers.size() < k) {
                updateMLinCenter();
                farthestPoint();
                updateMLinCenter();
                if (centers.size() < k) {
                    cannotLinkSet.sort((list1, list2) -> Integer.compare(list2.size(), list1.size()));
                    exchangeCenter();
                }
                updateMLinCenter();
            }

            if (centers.size() > k || (centers.size() == k && !judgeCannot(R))) {
                flag = false;
            }

            if (flag) {
                Rmax = (Rmin + Rmax) / 2.0f;
            } else {
                Rmin = (Rmin + Rmax) / 2.0f;
            }
        }
    }


    // Method to exchange centers
    private void exchangeCenter() {
        while (!judgeCannot(R)) {
            for (List<Integer> cannotLink : cannotLinkSet) {
                List<Integer> repeatLeft = new ArrayList<>();
                List<Integer> repeatRight = new ArrayList<>();
                boolean Flag_change = false;

                processCannotLinks(cannotLink, repeatLeft, repeatRight);
                centers.removeAll(repeatLeft);
                cannotLink.removeAll(repeatRight);

                if (!cannotLink.isEmpty()) {

                    Flag_change = matchinga.input(pointList, centers, cannotLink, mustLinkSet, R);
//                    Flag_change = RDS.input(pointList, centers, cannotLink, mustLinkSet, R);
                }
                centers.addAll(repeatLeft);
                cannotLink.addAll(repeatRight);
                updateMLinCenter();
                if (centers.size() > k || (Flag_change && judgeCannot(R))) {
                    break;
                }
            }
            if (centers.size() > k) {
                break;
            }
        }
    }

    private void processCannotLinks(List<Integer> cannotLink, List<Integer> repeatLeft, List<Integer> repeatRight) {

        for (int e : cannotLink) {
            Point point = pointList.get(e);
            int mustID = point.getMustID();
            if (centers.contains(e)) {
                repeatLeft.add(e);
                repeatRight.add(e);
            } else if (must.contains(mustID)) {
                repeatRight.add(e);
                for (Integer center : centers) {
                    if (pointList.get(center).getMustID() == mustID) {
                        repeatLeft.add(center);
                    }
                }
            }
        }
    }

    // Method to update must list
    private void updateMLinCenter() {
        must.clear();
        for (int i = 0; i < mustLinkSet.size(); i++) {
            ArrayList<Integer> commonElements = new ArrayList<>(mustLinkSet.get(i));
            commonElements.retainAll(centers);
            if (!commonElements.isEmpty()) {
                must.add(i);
            }
        }
    }

    // Method to update centers
    private void farthestPoint() {
        for (int i = 0; i < pointList.size(); i++) {
            Point currentPoint = pointList.get(i);
            int mustID = currentPoint.getMustID();
            boolean Flag_f = true;

            if (currentPoint.getConID() == -1 && !must.contains(mustID)) {
                List<Integer> linkedPoints = mustID != -1 ? mustLinkSet.get(mustID) : Collections.singletonList(i);
                int mark = i;
                for (int center : centers) {
                    float maxDist = 0;
                    for (int linkedPointID : linkedPoints) {
                        float dist = distance(pointList.get(center), pointList.get(linkedPointID));
                        if (maxDist < dist) {
                            maxDist = dist;
                            mark = linkedPointID;
                        }
                    }
                    if (maxDist <= R) {
                        Flag_f = false;
                        break;
                    }
                }

                if (Flag_f) {
                    centers.add(mark);
                    if (mustID != -1) must.add(mustID);
                }
            }
        }
    }

    private void getKPoint(List<Integer> input_center, Float[] Dist_Points) {
        int farPID = -1;
        float farPdistance = -1;

        for (int i = 0; i < pointList.size(); i++) {
            Point currentPoint = pointList.get(i);
            int mustID = currentPoint.getMustID();

            if (!must.contains(mustID)) {
                List<Integer> linkedPoints = mustID != -1 ? mustLinkSet.get(mustID) : Collections.singletonList(i);
                int mark = i;
                float minDist = Float.MAX_VALUE;

                for (int center : input_center) {
                    float maxDist = 0;
                    for (int linkedPointID : linkedPoints) {
                        float dist = distance(pointList.get(center), pointList.get(linkedPointID));
                        if (maxDist < dist) {
                            maxDist = dist;
                            mark = linkedPointID;
                        }
                    }
                    minDist = Math.min(minDist, maxDist);
                }
                for (int linkedPointID : linkedPoints) {
                    Dist_Points[linkedPointID] = Math.min(minDist, Dist_Points[linkedPointID]);
                }
                if (Dist_Points[mark] > farPdistance) {
                    farPdistance = minDist;
                    farPID = mark;
                }
            }
        }
        centers.add(farPID);
        if (pointList.get(farPID).getMustID() != -1) must.add(pointList.get(farPID).getMustID());
        if (centers.size() < k) getKPoint(Collections.singletonList(farPID), Dist_Points);
    }

    private String algorithm(Random r) {
        // Select the center processing
        long startTime1 = System.nanoTime();

        List<List<Integer>> candidateCenterSet = findCandidateCenterSet();
        List<Integer> initCenter = initializeCenter(candidateCenterSet, r);
        long startTime2 = System.nanoTime();

        performBinarySearch(initCenter);// Ensure R: binary search
        long time3 = System.nanoTime() - startTime2;

        if (centers.size() < k) {
            Float[] Dist_Points = new Float[N];
            Arrays.fill(Dist_Points, 0.0f);
            getKPoint(centers, Dist_Points);
        }

        setClusterIDs(centers);
        updateMLinCenter();

        assignGene();
        assignMust();
        assignCannot();

        return (System.nanoTime() - startTime1) + "," + time3;
    }


    private boolean judgeCannot(float inputR) {
        boolean flag = true;
        for (List<Integer> cannotLink : cannotLinkSet) {
            List<Integer> pointLeft = new ArrayList<>(centers);
            List<Integer> pointRight = new ArrayList<>(cannotLink);
            List<Integer> repeatLeft = new ArrayList<>();
            List<Integer> repeatRight = new ArrayList<>();
            processCannotLinks(cannotLink, repeatLeft, repeatRight);
            pointLeft.removeAll(repeatLeft);
            pointRight.removeAll(repeatRight);

            flag = judgeCannotOnce(inputR, pointLeft, pointRight);
            if (!flag) {
                return false;
            }
        }
        return flag;
    }


    private boolean judgeCannotOnce(float inputR, List<Integer> pointsLeft, List<Integer> pointsRight) {

        //for cannot-link:add center
        matching2 minMaxMatchingNew = new matching2(pointList, pointsLeft, pointsRight, inputR, mustLinkSet);
        return minMaxMatchingNew.searchcount() >= pointsRight.size();
    }


    private void assignMust() {

        for (int mustID : must) {
            for (Integer centerID : centers) {
                if (pointList.get(centerID).getMustID() == mustID) {
                    for (int linkedPointID : mustLinkSet.get(mustID)) {
                        pointList.get(linkedPointID).setClusterID(pointList.get(centerID).getClusterID());
                    }
                    break;
                }
            }
        }


        for (int j = 0; j < mustLinkSet.size(); j++) {
            int mark = -1;
            if (!must.contains(j)) {
                float minDist = Float.MAX_VALUE;

                for (int l = 0; l < centers.size(); l++) {
                    float maxDist = 0;
                    Point centerPoint = pointList.get(centers.get(l));

                    for (int linkedPointID : mustLinkSet.get(j)) {
                        float tempDist = Math.max(maxDist, distance(centerPoint, pointList.get(linkedPointID)));
                        if (centerPoint.getMustID() != -1) {
                            for (int mustID : mustLinkSet.get(centerPoint.getMustID())) {
                                tempDist = Math.max(tempDist, distance(pointList.get(linkedPointID), pointList.get(mustID)));
                            }
                        }
                        maxDist = Math.max(maxDist, tempDist);
                    }
                    if (minDist > maxDist) {
                        minDist = maxDist;
                        mark = l;
                    }
                    for (int p_i : mustLinkSet.get(j)) {
                        Point p = pointList.get(p_i);
                        p.setClusterID(mark);
                    }
                }
            }
        }
    }

    private void assignCannot() {

        //for cannot-link:add center
        for (List<Integer> cannotLink : cannotLinkSet) {
            ArrayList<Integer> pointsLeft = new ArrayList<>(centers);
            ArrayList<Integer> pointsRight = new ArrayList<>(cannotLink);
            List<Integer> repeatLeft = new ArrayList<>();
            List<Integer> repeatRight = new ArrayList<>();
            processCannotLinks(cannotLink, repeatLeft, repeatRight);
            pointsLeft.removeAll(repeatLeft);
            pointsRight.removeAll(repeatRight);

//            float inRmax = RMark[0];
            float inRmax = R;
            while (pointsRight.size() > 1) {

                float inRmin = 0;
                float inR = (inRmin + inRmax) / 2.0f;
                boolean conditionMet = false;
                while (inRmax - inR > errorBar || !conditionMet) {
                    inR = (inRmin + inRmax) / 2.0f;
                    if (inR == inRmax || inR == inRmin) {
                        inR = inRmax;
                        inRmin = inRmax;
//                        break;
                    }
                    if (judgeCannotOnce(inR, pointsLeft, pointsRight)) {
                        inRmax = inR;
                        conditionMet = true;
                    } else {
                        inRmin = inR;
                        conditionMet = false;
                    }
                }
                matching2 MAlg = new matching2(pointList, pointsLeft, pointsRight, inR, mustLinkSet);
                MAlg.assignCannot(pointList, mustLinkSet);

                int markFar = 0;
                float farPointDist = 0;
                for (int j = 0; j < pointsRight.size(); j++) {
                    float dist = distance(pointList.get(pointsRight.get(j)), pointList.get(pointsLeft.get(MAlg.match[j])));
                    if (dist > farPointDist) {
                        markFar = j;
                        farPointDist = dist;
                    }
                }
                pointsLeft.remove(MAlg.match[markFar]);
                pointsRight.remove(markFar);
            }

            if (pointsRight.size() == 1) {
                Point rightPoint = pointList.get(pointsRight.get(0));
                int mustID = rightPoint.getMustID();

                float minDist = Float.MAX_VALUE;
                int mark = 0;

                for (int j = 0; j < pointsLeft.size(); j++) {
                    Point leftPoint = pointList.get(pointsLeft.get(j));
                    float dmax = distance(rightPoint, leftPoint);

                    if (mustID != -1) {
                        for (int linkedPointID : mustLinkSet.get(mustID)) {
                            Point linkedPoint = pointList.get(linkedPointID);
                            dmax = Math.max(dmax, distance(linkedPoint, leftPoint));
                        }
                    }
                    if (dmax < minDist) {
                        mark = j;
                        minDist = dmax;
                    }
                }
                int clusterID = pointList.get(pointsLeft.get(mark)).getClusterID();
                rightPoint.setClusterID(clusterID);
                if (mustID != -1) {
                    for (int linkedPointID : mustLinkSet.get(mustID)) {
                        pointList.get(linkedPointID).setClusterID(clusterID);
                    }
                }
            }
        }
    }

    private static void assignGene() {
        for (Point point : pointList) {
            if (point.getMustID() == -1 && point.getConID() == -1) {
                float mindistance = Float.MAX_VALUE;
                for (int l = 0; l < centers.size(); l++) {
                    float dist = tools.distance(pointList.get(centers.get(l)), point);
                    if (dist < mindistance) {
                        point.setClusterID(l);
                        mindistance = dist;
                    }
                }
            }
        }
    }
}





//
//
//
//
//
//package OurAlg_impro_matching;
//
//import Tools.*;
//
//import java.util.*;
//
//import static Tools.tools.distance;
//
//public class Approx_improv_matching_kCenter {
//
//    public static int k, N;
//    public static ArrayList<Point> pointList;
//    public static ArrayList<ArrayList<Integer>> cannotLinkSet, mustLinkSet;
//    static float[] RMark;
//    public static long time3 = 0;
//    static float R, errorBar, Rmax, Rmin;
//    static ArrayList<Integer> must = new ArrayList<>(), centers;
//
//    public Approx_improv_matching_kCenter(ArrayList<Point> pointList, ArrayList<ArrayList<Integer>> CannotList, ArrayList<ArrayList<Integer>> MustList, int k, float[] RMark) {
//        this.pointList = pointList;
//        cannotLinkSet = CannotList;
//        mustLinkSet = MustList;
//        this.k = k;
//        this.RMark = RMark;
//    }
//
//    // Method to set clusterID for centers
//    private void setClusterIDs(List<Integer> centers) {
//        for (int i = 0; i < centers.size(); i++) {
//            pointList.get(centers.get(i)).setClusterID(i);
//        }
//    }
//
//    public float[] vertexCover(Random r) {
//        System.out.println("Approx_improv_matching_kCenter:___________________________________________");
//        N = pointList.size();
//        float[] accuracy;
//        time3 = 0;
//        Rmax = RMark[0];
//        Rmin = RMark[1];
//        errorBar = RMark[2];
//        centers = new ArrayList<>();
//        must.clear();
//
//        for (int i = 0; i < N; i++) {
//            pointList.get(i).setClusterID(-1);
//        }
//
//        String times = algorithm(r);
//        accuracy = Metrics.metrics(N, pointList, centers, k, times + ", " + time3);
//        return accuracy;
//    }
//
//    // Method to find the candidate center set
//    private List<List<Integer>> findCandidateCenterSet() {
//        List<List<Integer>> candidateCenterSet = new ArrayList<>();
//        int maximumSize = 0;
//        for (ArrayList<Integer> integerArrayList : cannotLinkSet) {
//            maximumSize = Math.max(maximumSize, integerArrayList.size());
//        }
//        for (ArrayList<Integer> integers : cannotLinkSet) {
//            if (integers.size() == maximumSize) {
//                candidateCenterSet.add(integers);
//            }
//        }
//        return candidateCenterSet;
//    }
//
//    // Method to initialize center
//    private List<Integer> initializeCenter(List<List<Integer>> candidateCenterSet, Random r) {
//        List<Integer> initCenter = new ArrayList<>();
//        if (candidateCenterSet.size() == 0) {
//            initCenter.add(r.nextInt(N));
//        } else {
//            initCenter.addAll(candidateCenterSet.get(r.nextInt(candidateCenterSet.size())));
//        }
//        return initCenter;
//    }
//
//    // Method to perform binary search to ensure R
//    private void performBinarySearch(List<Integer> initCenter) {
//
//        boolean flag = false;
//        while (Rmax - Rmin > errorBar || !flag) {
//            R = (Rmin + Rmax) / 2.0f;
//            if (R == Rmin || R == Rmax) {
//                Rmin = R = Rmax;
//                break;
//            }
//            flag = true;
//            centers.clear();
//            centers.addAll(initCenter);
//
//            // Exchange the center
//            if (centers.size() < k) {
//                System.out.println("?");
//                updateMLinCenter();
//                farthestPoint();
////                if (centers.size() < k) {
//                    cannotLinkSet.sort((list1, list2) -> Integer.compare(list2.size(), list1.size()));
//                    exchangeCenter();
////                }
//                updateMLinCenter();
//            }
//
//            if (centers.size() > k || (centers.size() == k && !judgeCannot(R))) {
//                flag = false;
//            }
//
//            if (flag) {
//                Rmax = (Rmin + Rmax) / 2.0f;
//            } else {
//                Rmin = (Rmin + Rmax) / 2.0f;
//            }
//        }
//    }
//
//
//    // Method to exchange centers
//    private void exchangeCenter() {
//        while (!judgeCannot(R)) {
//            for (List<Integer> cannotLink : cannotLinkSet) {
//                List<Integer> repeatLeft = new ArrayList<>();
//                List<Integer> repeatRight = new ArrayList<>();
//                boolean Flag_change = false;
//
//                processCannotLinks(cannotLink, repeatLeft, repeatRight);
//                centers.removeAll(repeatLeft);
//                cannotLink.removeAll(repeatRight);
//                if (!cannotLink.isEmpty()) {
//                    Flag_change = matchinga.input(pointList, centers, cannotLink, mustLinkSet, R);
//                }
//                centers.addAll(repeatLeft);
//                cannotLink.addAll(repeatRight);
//                updateMLinCenter();
//                if (centers.size() > k || (Flag_change && judgeCannot(R))) {
//                    break;
//                }
//            }
//
//            if (centers.size() > k) {
//                break;
//            }
//        }
//    }
//
//    private void processCannotLinks(List<Integer> cannotLink, List<Integer> repeatLeft, List<Integer> repeatRight) {
//
//        for (int e : cannotLink) {
//            Point point = pointList.get(e);
//            int mustID = point.getMustID();
//            if (centers.contains(e)) {
//                repeatLeft.add(e);
//                repeatRight.add(e);
//            } else if (must.contains(mustID)) {
//                repeatRight.add(e);
//                for (Integer center : centers) {
//                    if (pointList.get(center).getMustID() == mustID) {
//                        repeatLeft.add(center);
//                    }
//                }
//            }
//        }
//    }
//
//    // Method to update must list
//    private void updateMLinCenter() {
//        must.clear();
//        for (int i = 0; i < mustLinkSet.size(); i++) {
//            ArrayList<Integer> commonElements = new ArrayList<>(mustLinkSet.get(i));
//            commonElements.retainAll(centers);
//            if (!commonElements.isEmpty()) {
//                must.add(i);
//            }
//        }
//    }
//
//    // Method to update centers
//    private void farthestPoint() {
//        for (int i = 0; i < pointList.size(); i++) {
//            Point currentPoint = pointList.get(i);
//            int mustID = currentPoint.getMustID();
//            boolean Flag_f = true;
//
//            if (currentPoint.getConID() == -1 && !must.contains(mustID)) {
//                List<Integer> linkedPoints = mustID != -1 ? mustLinkSet.get(mustID) : Collections.singletonList(i);
//                int mark = i;
//                for (int center : centers) {
//                    float maxDist = 0;
//                    for (int linkedPointID : linkedPoints) {
//                        float dist = distance(pointList.get(center), pointList.get(linkedPointID));
//                        if (maxDist < dist) {
//                            maxDist = dist;
//                            mark = linkedPointID;
//                        }
//                    }
//                    if (maxDist <= R) {
//                        Flag_f = false;
//                        break;
//                    }
//                }
//
//                if (Flag_f) {
//                    centers.add(mark);
//                    if (mustID != -1) must.add(mustID);
//                }else {
//                    break;
//                }
//            }
//        }
//    }
//
//    private void getKPoint(List<Integer> input_center, Float[] Dist_Points) {
//        int farPID = -1;
//        float farPdistance = -1;
//
//        for (int i = 0; i < pointList.size(); i++) {
//            Point currentPoint = pointList.get(i);
//            int mustID = currentPoint.getMustID();
//
//            if (!must.contains(mustID)) {
//                List<Integer> linkedPoints = mustID != -1 ? mustLinkSet.get(mustID) : Collections.singletonList(i);
//                int mark = i;
//                float minDist = Float.MAX_VALUE;
//
//                for (int center : input_center) {
//                    float maxDist = 0;
//                    for (int linkedPointID : linkedPoints) {
//                        float dist = distance(pointList.get(center), pointList.get(linkedPointID));
//                        if (maxDist < dist) {
//                            maxDist = dist;
//                            mark = linkedPointID;
//                        }
//                    }
//                    minDist = Math.min(minDist, maxDist);
//                }
//                for (int linkedPointID : linkedPoints) {
//                    Dist_Points[linkedPointID] = Math.min(minDist, Dist_Points[linkedPointID]);
//                }
//                if (Dist_Points[mark] > farPdistance) {
//                    farPdistance = minDist;
//                    farPID = mark;
//                }
//            }
//        }
//        centers.add(farPID);
//        if (pointList.get(farPID).getMustID() != -1) must.add(pointList.get(farPID).getMustID());
//        if (centers.size() < k) getKPoint(Collections.singletonList(farPID), Dist_Points);
//    }
//
//    private String algorithm(Random r) {
//        // Select the center processing
//        long startTime1 = System.nanoTime();
//
//        List<List<Integer>> candidateCenterSet = findCandidateCenterSet();
//        List<Integer> initCenter = initializeCenter(candidateCenterSet, r);
//        long startTime2 = System.nanoTime();
//        performBinarySearch(initCenter);// Ensure R: binary search
//        long time3 = System.nanoTime() - startTime2;
//
//        if (centers.size() < k) {
//            Float[] Dist_Points = new Float[N];
//            Arrays.fill(Dist_Points, 0.0f);
//            getKPoint(centers, Dist_Points);
//        }
//
//        setClusterIDs(centers);
//        updateMLinCenter();
//
//        assignGene();
//        assignMust();
//        assignCannot();
//
//        return (System.nanoTime() - startTime1) + "," + time3;
//    }
//
//
//    private boolean judgeCannot(float inputR) {
//        boolean flag = true;
//        for (List<Integer> cannotLink : cannotLinkSet) {
//            List<Integer> pointLeft = new ArrayList<>(centers);
//            List<Integer> pointRight = new ArrayList<>(cannotLink);
//            List<Integer> repeatLeft = new ArrayList<>();
//            List<Integer> repeatRight = new ArrayList<>();
//            processCannotLinks(cannotLink, repeatLeft, repeatRight);
//            pointLeft.removeAll(repeatLeft);
//            pointRight.removeAll(repeatRight);
//
//            flag = judgeCannotOnce(inputR, pointLeft, pointRight);
//            if (!flag) {
//                return false;
//            }
//        }
//        return flag;
//    }
//
//
//    private boolean judgeCannotOnce(float inputR, List<Integer> pointsLeft, List<Integer> pointsRight) {
//
//        //for cannot-link:add center
//        matching2 minMaxMatchingNew = new matching2(pointList, pointsLeft, pointsRight, inputR, mustLinkSet);
//        return minMaxMatchingNew.searchcount() >= pointsRight.size();
//    }
//
//
//    private void assignMust() {
//
//        for (int mustID : must) {
//            for (Integer centerID : centers) {
//                if (pointList.get(centerID).getMustID() == mustID) {
//                    for (int linkedPointID : mustLinkSet.get(mustID)) {
//                        pointList.get(linkedPointID).setClusterID(pointList.get(centerID).getClusterID());
//                    }
//                    break;
//                }
//            }
//        }
//
//
//        for (int j = 0; j < mustLinkSet.size(); j++) {
//            int mark = -1;
//            if (!must.contains(j)) {
//                float minDist = Float.MAX_VALUE;
//
//                for (int l = 0; l < centers.size(); l++) {
//                    float maxDist = 0;
//                    Point centerPoint = pointList.get(centers.get(l));
//
//                    for (int linkedPointID : mustLinkSet.get(j)) {
//                        float tempDist = Math.max(maxDist, distance(centerPoint, pointList.get(linkedPointID)));
//                        if (centerPoint.getMustID() != -1) {
//                            for (int mustID : mustLinkSet.get(centerPoint.getMustID())) {
//                                tempDist = Math.max(tempDist, distance(pointList.get(linkedPointID), pointList.get(mustID)));
//                            }
//                        }
//                        maxDist = Math.max(maxDist, tempDist);
//                    }
//                    if (minDist > maxDist) {
//                        minDist = maxDist;
//                        mark = l;
//                    }
//                    for (int p_i : mustLinkSet.get(j)) {
//                        Point p = pointList.get(p_i);
//                        p.setClusterID(mark);
//                    }
//                }
//            }
//        }
//    }
//
//    private void assignCannot() {
//
//        //for cannot-link:add center
//        for (List<Integer> cannotLink : cannotLinkSet) {
//            ArrayList<Integer> pointsLeft = new ArrayList<>(centers);
//            ArrayList<Integer> pointsRight = new ArrayList<>(cannotLink);
//            List<Integer> repeatLeft = new ArrayList<>();
//            List<Integer> repeatRight = new ArrayList<>();
//            processCannotLinks(cannotLink, repeatLeft, repeatRight);
//            pointsLeft.removeAll(repeatLeft);
//            pointsRight.removeAll(repeatRight);
//
//            float inRmax = R;
//            while (pointsRight.size() > 1) {
//
//                float inRmin = 0;
//                float inR = (inRmin + inRmax) / 2.0f;
//                boolean conditionMet = false;
//                while (inRmax - inR > errorBar || !conditionMet) {
//                    inR = (inRmin + inRmax) / 2.0f;
//                    if (inR == inRmax || inR == inRmin) {
//                        inR = inRmax;
////                        inRmin = inRmax;
//                        break;
//                    }
//                    if (judgeCannotOnce(inR, pointsLeft, pointsRight)) {
//                        inRmax = inR;
//                        conditionMet = true;
//                    } else {
//                        inRmin = inR;
//                        conditionMet = false;
//                    }
//                }
//                matching2 MAlg = new matching2(pointList, pointsLeft, pointsRight, inR, mustLinkSet);
//                MAlg.assignCannot(pointList, mustLinkSet);
//
//                int markFar = 0;
//                float farPointDist = 0;
//                for (int j = 0; j < pointsRight.size(); j++) {
//                    float dist = distance(pointList.get(pointsRight.get(j)), pointList.get(pointsLeft.get(MAlg.match[j])));
//                    if (dist > farPointDist) {
//                        markFar = j;
//                        farPointDist = dist;
//                    }
//                }
//                pointsLeft.remove(MAlg.match[markFar]);
//                pointsRight.remove(markFar);
//            }
//
//            if (pointsRight.size() == 1) {
//                Point rightPoint = pointList.get(pointsRight.get(0));
//                int mustID = rightPoint.getMustID();
//
//                float minDist = Float.MAX_VALUE;
//                int mark = 0;
//
//                for (int j = 0; j < pointsLeft.size(); j++) {
//                    Point leftPoint = pointList.get(pointsLeft.get(j));
//                    float dmax = distance(rightPoint, leftPoint);
//
//                    if (mustID != -1) {
//                        for (int linkedPointID : mustLinkSet.get(mustID)) {
//                            Point linkedPoint = pointList.get(linkedPointID);
//                            dmax = Math.max(dmax, distance(linkedPoint, leftPoint));
//                        }
//                    }
//                    if (dmax < minDist) {
//                        mark = j;
//                        minDist = dmax;
//                    }
//                }
//                int clusterID = pointList.get(pointsLeft.get(mark)).getClusterID();
//                rightPoint.setClusterID(clusterID);
//                if (mustID != -1) {
//                    for (int linkedPointID : mustLinkSet.get(mustID)) {
//                        pointList.get(linkedPointID).setClusterID(clusterID);
//                    }
//                }
//            }
//        }
//    }
//
//    private static void assignGene() {
//        for (Point point : pointList) {
//            if (point.getMustID() == -1 && point.getConID() == -1) {
//                float mindistance = Float.MAX_VALUE;
//                for (int l = 0; l < centers.size(); l++) {
//                    float dist = tools.distance(pointList.get(centers.get(l)), point);
//                    if (dist < mindistance) {
//                        point.setClusterID(l);
//                        mindistance = dist;
//                    }
//                }
//            }
//        }
//    }
//}
//
//
////package OurAlg_impro_matching;
////
////import Tools.*;
////
////import java.util.*;
////
////import static Tools.tools.distance;
////
////
////public class Approx_improv_matching_kCenter {
////
////    public static int k, N;
////    public static long time3;
////    public static boolean break_condi;
////    public static ArrayList<Point> pointList;
////    public static ArrayList<ArrayList<Integer>> cannotLinkSet, mustLinkSet;
////    static float[] RMark;
////    static float R, errorBar, Rmax, Rmin;
////    static ArrayList<Integer> must = new ArrayList<>(), centers;
////
////    public Approx_improv_matching_kCenter(ArrayList<Point> pointList, ArrayList<ArrayList<Integer>> CannotList, ArrayList<ArrayList<Integer>> MustList, int k, float[] RMark) {
////        this.pointList = pointList;
////        cannotLinkSet = CannotList;
////        mustLinkSet = MustList;
////        this.k = k;
////        this.RMark = RMark;
////    }
////
////    // Method to set clusterID for centers
////    private void setClusterIDs(List<Integer> centers) {
////        for (int i = 0; i < centers.size(); i++) {
////            pointList.get(centers.get(i)).setClusterID(i);
////        }
////    }
////
////    public float[] vertexCover(Random r) {
////        System.out.println("Approx_improv_matching_kCenter:___________________________________________");
////        N = pointList.size();
////        float[] accuracy;
////
////        Rmax = RMark[0];
////        Rmin = RMark[1];
////        errorBar = RMark[2];
////        centers = new ArrayList<>();
////        must.clear();
////
////        for (int i = 0; i < N; i++) {
////            pointList.get(i).setClusterID(-1);
////        }
////
////        String times = algorithm(r);
////        accuracy = Metrics.metrics(N, pointList, centers, k, times);
////
////        return accuracy;
////    }
////
////    // Method to find the candidate center set
////    private List<List<Integer>> findCandidateCenterSet() {
////        List<List<Integer>> candidateCenterSet = new ArrayList<>();
////        int maximumSize = 0;
////        for (ArrayList<Integer> integerArrayList : cannotLinkSet) {
////            maximumSize = Math.max(maximumSize, integerArrayList.size());
////        }
////        for (ArrayList<Integer> integers : cannotLinkSet) {
////            if (integers.size() == maximumSize) {
////                candidateCenterSet.add(integers);
////            }
////        }
////        return candidateCenterSet;
////    }
////
////    // Method to initialize center
////    private List<Integer> initializeCenter(List<List<Integer>> candidateCenterSet, Random r) {
////        List<Integer> initCenter = new ArrayList<>();
////        if (candidateCenterSet.size() == 0) {
////            initCenter.add(r.nextInt(N));
////        } else {
////            initCenter.addAll(candidateCenterSet.get(r.nextInt(candidateCenterSet.size())));
////        }
////        return initCenter;
////    }
////
////    // Method to perform binary search to ensure R
////    private void performBinarySearch(List<Integer> initCenter) {
////        boolean flag = false;
////        while (Rmax - Rmin > errorBar || !flag) {
////            R = (Rmin + Rmax) / 2.0f;
////            if (R == Rmin || R == Rmax) {
////                Rmin = R = Rmax;
////            }
////            flag = true;
////            centers.clear();
////            centers.addAll(initCenter);
////            updateMLinCenter();
////
////            farthestPoint();
////            // Exchange the center
////
////            if (centers.size() < k) {
////                cannotLinkSet.sort((list1, list2) -> Integer.compare(list2.size(), list1.size()));
////                exchangeCenter();
////                updateMLinCenter();
////            }
////
////
////            if (centers.size() > k || (centers.size() == k && !judgeCannot(R))) {
////                flag = false;
////            }
////
////            if (flag) {
////                Rmax = (Rmin + Rmax) / 2.0f;
////            } else {
////                Rmin = (Rmin + Rmax) / 2.0f;
////            }
////
////        }
////
////    }
////
////
////    // Method to exchange centers
////    private void exchangeCenter() {
//////        outer:
////        while (!judgeCannot(R)) {
////            for (List<Integer> cannotLink : cannotLinkSet) {
////                List<Integer> repeatLeft = new ArrayList<>();
////                List<Integer> repeatRight = new ArrayList<>();
//////                matchinga ma = new matchinga();
////                boolean Flag_change = false;
////
////                updateMLinCenter();
////                processCannotLinks(cannotLink, repeatLeft, repeatRight);
////                centers.removeAll(repeatLeft);
////                cannotLink.removeAll(repeatRight);
////
////                if (!cannotLink.isEmpty()) {
////                    Flag_change =  matchinga.input(pointList, centers, cannotLink, mustLinkSet, R,time3);
//////                    updateMLinCenter();
////                }
////
////                centers.addAll(repeatLeft);
////                cannotLink.addAll(repeatRight);
////                if ( centers.size() > k || Flag_change && judgeCannot(R)){
////                    break;
////                }
////            }
////            if (centers.size() > k) {
////                break;
////            }
////        }
////    }
////
////    private void processCannotLinks(List<Integer> cannotLink, List<Integer> repeatLeft, List<Integer> repeatRight) {
////
////        for (int e : cannotLink) {
////            Point point = pointList.get(e);
////            int mustID = point.getMustID();
////            if (centers.contains(e)) {
////                repeatLeft.add(e);
////                repeatRight.add(e);
////            } else if (must.contains(mustID)) {
////                repeatRight.add(e);
////                for (Integer center : centers) {
////                    if (pointList.get(center).getMustID() == mustID) {
////                        repeatLeft.add(center);
////                    }
////                }
////            }
////        }
////    }
////
////    // Method to update must list
////    private void updateMLinCenter() {
////        must.clear();
////        for (int i = 0; i < mustLinkSet.size(); i++) {
////            ArrayList<Integer> commonElements = new ArrayList<>(mustLinkSet.get(i));
////            commonElements.retainAll(centers);
////            if (!commonElements.isEmpty()) {
////                must.add(i);
////            }
////        }
////    }
////
////    // Method to update centers
////    private void farthestPoint() {
////        for (int i = 0; i < pointList.size(); i++) {
////            Point currentPoint = pointList.get(i);
////            int mustID = currentPoint.getMustID();
////            boolean Flag_f = true;
////
////            if (currentPoint.getConID() == -1 && !must.contains(mustID)) {
////                List<Integer> linkedPoints = mustID != -1 ? mustLinkSet.get(mustID) : Collections.singletonList(i);
////                int mark = i;
////                for (int center : centers) {
////                    float maxDist = 0;
////                    for (int linkedPointID : linkedPoints) {
////                        float dist = distance(pointList.get(center), pointList.get(linkedPointID));
////                        if (maxDist < dist) {
////                            maxDist = dist;
////                            mark = linkedPointID;
////                        }
////                    }
////                    if (maxDist <= R) {
////                        Flag_f = false;
////                        break;
////                    }
////                }
////
////                if (Flag_f) {
////                    centers.add(mark);
////                    if (mustID != -1) must.add(mustID);
////                }
////            }
////        }
////    }
////
////    private void getKPoint(List<Integer> input_center, Float[] Dist_Points) {
////        int farPID = -1;
////        float farPdistance = -1;
////
////        for (int i = 0; i < pointList.size(); i++) {
////            Point currentPoint = pointList.get(i);
////            int mustID = currentPoint.getMustID();
////
////            if (!must.contains(mustID)) {
////                List<Integer> linkedPoints = mustID != -1 ? mustLinkSet.get(mustID) : Collections.singletonList(i);
////                int mark = i;
////                float minDist = Float.MAX_VALUE;
////
////                for (int center : input_center) {
////                    float maxDist = 0;
//////                    int markMaxDist = -1;
////                    for (int linkedPointID : linkedPoints) {
////                        float dist = distance(pointList.get(center), pointList.get(linkedPointID));
////                        if (maxDist < dist) {
////                            maxDist = dist;
////                            mark = linkedPointID;
////                        }
////                    }
////                    minDist = Math.min(minDist, maxDist);
////                }
////                for (int linkedPointID : linkedPoints) {
////                    Dist_Points[linkedPointID] = Math.min(minDist, Dist_Points[linkedPointID]);
////                }
////                if (Dist_Points[mark] > farPdistance) {
////                    farPdistance = minDist;
////                    farPID = mark;
////                }
////            }
////        }
////        centers.add(farPID);
////        if (pointList.get(farPID).getMustID() != -1) must.add(pointList.get(farPID).getMustID());
////        if (centers.size() < k) getKPoint(Collections.singletonList(farPID), Dist_Points);
////    }
////
////    private String algorithm(Random r) {
////        long startTime1 = System.nanoTime();
////        // Find the center
////        List<List<Integer>> candidateCenterSet = findCandidateCenterSet();
////        List<Integer> initCenter = initializeCenter(candidateCenterSet, r);
////
////        performBinarySearch(initCenter);// Ensure R: binary search
////        long time3 = System.nanoTime() - startTime1;
////
////        if (centers.size() < k) {
////            Float[] Dist_Points = new Float[N];
////            Arrays.fill(Dist_Points, 0.0f);
////            getKPoint(centers, Dist_Points);
////        }
////
////        setClusterIDs(centers);
////        updateMLinCenter();
////        assignGene();
////        assignMust();
////        assignCannot();
////
////        return (System.nanoTime() - startTime1) + "," + time3;
////    }
////
////
////    private boolean judgeCannot(float inputR) {
////        boolean flag = true;
////        for (List<Integer> cannotLink : cannotLinkSet) {
////            List<Integer> pointLeft = new ArrayList<>(centers);
////            List<Integer> pointRight = new ArrayList<>(cannotLink);
////            List<Integer> repeatLeft = new ArrayList<>();
////            List<Integer> repeatRight = new ArrayList<>();
////            processCannotLinks(cannotLink, repeatLeft, repeatRight);
////            pointLeft.removeAll(repeatLeft);
////            pointRight.removeAll(repeatRight);
////
////            flag = judgeCannotOnce(inputR, pointLeft, pointRight);
////            if (!flag) {
////                return false;
////            }
////        }
////        return flag;
////    }
////
////
////    private boolean judgeCannotOnce(float inputR, List<Integer> pointsLeft, List<Integer> pointsRight) {
////
////        //for cannot-link:add center
////        matching2 minMaxMatchingNew = new matching2(pointList, pointsLeft, pointsRight, inputR, mustLinkSet);
////        return minMaxMatchingNew.searchcount() >= pointsRight.size();
////    }
////
////
////    private void assignMust() {
////
////        for (int mustID : must) {
////            for (Integer centerID : centers) {
////                if (pointList.get(centerID).getMustID() == mustID) {
////                    for (int linkedPointID : mustLinkSet.get(mustID)) {
////                        pointList.get(linkedPointID).setClusterID(pointList.get(centerID).getClusterID());
////                    }
////                    break;
////                }
////            }
////        }
////
////
////        for (int j = 0; j < mustLinkSet.size(); j++) {
////            int mark = -1;
////            if (!must.contains(j)) {
////                float minDist = Float.MAX_VALUE;
////
////                for (int l = 0; l < centers.size(); l++) {
////                    float maxDist = 0;
////                    Point centerPoint = pointList.get(centers.get(l));
////
////                    for (int linkedPointID : mustLinkSet.get(j)) {
////                        float tempDist = Math.max(maxDist, distance(centerPoint, pointList.get(linkedPointID)));
////                        if (centerPoint.getMustID() != -1) {
////                            for (int mustID : mustLinkSet.get(centerPoint.getMustID())) {
////                                tempDist = Math.max(tempDist, distance(pointList.get(linkedPointID), pointList.get(mustID)));
////                            }
////                        }
////                        maxDist = Math.max(maxDist, tempDist);
////                    }
////                    if (minDist > maxDist) {
////                        minDist = maxDist;
////                        mark = l;
////                    }
////                    for (int p_i : mustLinkSet.get(j)) {
////                        Point p = pointList.get(p_i);
////                        p.setClusterID(mark);
////                    }
////                }
////            }
////        }
////    }
////
////    private void assignCannot() {
////
////        //for cannot-link:add center
////        for (List<Integer> cannotLink : cannotLinkSet) {
////            ArrayList<Integer> pointsLeft = new ArrayList<>(centers);
////            ArrayList<Integer> pointsRight = new ArrayList<>(cannotLink);
////            List<Integer> repeatLeft = new ArrayList<>();
////            List<Integer> repeatRight = new ArrayList<>();
////            processCannotLinks(cannotLink, repeatLeft, repeatRight);
////            pointsLeft.removeAll(repeatLeft);
////            pointsRight.removeAll(repeatRight);
////
////            float inRmax = R;
////            while (pointsRight.size() > 1) {
////
////                float inRmin = 0;
////                float inR = R;
////                boolean conditionMet = false;
////                while (inRmax - inR > errorBar || !conditionMet) {
////
////                    inR = (inRmin + inRmax) / 2.0f;
////                    if (inR == inRmax || inR == inRmin) {
////                        inR = inRmax;
////                        inRmin = inRmax;
////                    }
////                    if (judgeCannotOnce(inR, pointsLeft, pointsRight)) {
////                        inRmax = inR;
////                        conditionMet = true;
////                    } else {
////                        inRmin = inR;
////                        conditionMet = false;
////                    }
////                }
////                matching2 MAlg = new matching2(pointList, pointsLeft, pointsRight, inR, mustLinkSet);
////                MAlg.assignCannot(pointList, mustLinkSet);
////
////                int markFar = 0;
////                float farPointDist = 0;
////                for (int j = 0; j < pointsRight.size(); j++) {
////                    float dist = distance(pointList.get(pointsRight.get(j)), pointList.get(pointsLeft.get(MAlg.match[j])));
////                    if (dist > farPointDist) {
////                        markFar = j;
////                        farPointDist = dist;
////                    }
////                }
////                pointsLeft.remove(MAlg.match[markFar]);
////                pointsRight.remove(markFar);
////            }
////
////            if (pointsRight.size() == 1) {
////                Point rightPoint = pointList.get(pointsRight.get(0));
////                int mustID = rightPoint.getMustID();
////
////                float minDist = Float.MAX_VALUE;
////                int mark = 0;
////
////                for (int j = 0; j < pointsLeft.size(); j++) {
////                    Point leftPoint = pointList.get(pointsLeft.get(j));
////                    float dmax = distance(rightPoint, leftPoint);
////
////                    if (mustID != -1) {
////                        for (int linkedPointID : mustLinkSet.get(mustID)) {
////                            Point linkedPoint = pointList.get(linkedPointID);
////                            dmax = Math.max(dmax, distance(linkedPoint, leftPoint));
////                        }
////                    }
////                    if (dmax < minDist) {
////                        mark = j;
////                        minDist = dmax;
////                    }
////                }
////                int clusterID = pointList.get(pointsLeft.get(mark)).getClusterID();
////                rightPoint.setClusterID(clusterID);
////                if (mustID != -1) {
////                    for (int linkedPointID : mustLinkSet.get(mustID)) {
////                        pointList.get(linkedPointID).setClusterID(clusterID);
////                    }
////                }
////            }
////        }
////    }
////
////    private static void assignGene() {
////        for (Point point : pointList) {
////            if (point.getMustID() == -1 && point.getConID() == -1) {
////                float mindistance = Float.MAX_VALUE;
////                for (int l = 0; l < centers.size(); l++) {
////                    float dist = tools.distance(pointList.get(centers.get(l)), point);
////                    if (dist < mindistance) {
////                        point.setClusterID(l);
////                        mindistance = dist;
////                    }
////                }
////            }
////        }
////    }
////
////}
//////
//////
//////    private int farthestPoint(boolean addCenter) {
//////        int farPID = -1;
//////        float farPdistance = -1;
//////
//////        for (int i = 0; i < pointList.size(); i++) {
//////            Point currentPoint = pointList.get(i);
//////            int mustID = currentPoint.getMustID();
//////
//////            if ((!addCenter || currentPoint.getConID() == -1) && !must.contains(mustID)) {
//////                List<Integer> linkedPoints = mustID != -1 ? mustLinkSet.get(mustID) : Collections.emptyList();
//////                int mark = i;
//////                float minDist = Float.MAX_VALUE;
//////
//////                if (mustID != -1) {
//////                    for (int center : centers) {
//////                        float maxDist = 0;
//////                        int markMaxDist = -1;
//////                        for (int linkedPointID : linkedPoints) {
//////                            float dist = distance(pointList.get(center), pointList.get(linkedPointID));
//////                            if (maxDist < dist) {
//////                                maxDist = dist;
//////                                markMaxDist = linkedPointID;
//////                            }
//////                        }
//////                        if (minDist > maxDist) {
//////                            minDist = maxDist;
//////                            mark = markMaxDist;
//////                        }
//////                    }
//////                } else {
//////                    for (int center : centers) {
//////                        float maxDist = distance(currentPoint, pointList.get(center));
//////                        minDist = Math.min(maxDist, minDist);
//////                    }
//////                }
//////
//////                if (!addCenter && minDist > farPdistance) {
//////                    farPdistance = minDist;
//////                    farPID = mark;
//////                } else if (addCenter && minDist > R) {
//////                    centers.add(mark);
//////                    if (mustID != -1) {
//////                        must.add(mustID);
//////                    }
//////                }
//////            }
//////        }
//////        return farPID;
//////    }
//
//
