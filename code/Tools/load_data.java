package Tools;

import addConstraints.OutPut;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class load_data {
    public static ArrayList<Point> readMyFile(String filePath, int index) {

        if (filePath.contains("wine") || filePath.contains("cnae")) {
            filePath = "/Datasets/" + filePath + ".data"; 
            int markPosition = 0;
            return loadWineData(filePath, markPosition, "[,]");
        } else if (filePath.contains("skin")) {
            filePath = "/Datasets/" + filePath + ".data"; 
            int markPosition = 3;
            return loadWineData(filePath, markPosition, "[\\s+]");
        } else if (filePath.contains("covertype")) {
            filePath =  "/Datasets/" + filePath + ".data"; 
            int markPosition = 54;
            return loadWineData(filePath, markPosition, "[,]");
        } else if (filePath.contains("wide09")) {
            filePath = "/Datasets/" + filePath + ".data"; 
            int markPosition = 21;
            return loadWineData(filePath, markPosition, "[,]");
        } else if (filePath.contains("kdd")) {
            filePath =  "/Datasets/" + filePath + ".data"; 
            int markPosition = 41;
            return loadWineData(filePath, markPosition, "[,]");
        }else if (filePath.contains("simu")) {
            filePath =  "/Datasets/" + filePath + ".data"; 
            int markPosition = 50;
            return loadWineData(filePath, markPosition, "[,]");
        }
        return null;
    }

    public static ArrayList<Point> loadWineData(String filePath, int markPosition, String delimiter) {
        ArrayList<Point> pointList = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            int i = 0;
            Set<Integer> label_numbers = new HashSet<>();
//            Map<String, Integer> countMap = new HashMap<>();
            while ((line = br.readLine()) != null) {
                String[] values = Arrays.stream(line.split(delimiter))
                        .map(String::trim)
                        .filter(s -> !s.isEmpty())
                        .toArray(String[]::new);
                int mark = (int) Double.parseDouble(values[markPosition].replaceAll("[\\[\\] ]", "").trim());
                int d = values.length - 1;
                float[] positionPoint = new float[d];

                for (int j = 0; j < d; j++) {
                    positionPoint[j] = Float.parseFloat(values[j + (markPosition + 1) % (d + 1)].replaceAll("[\\[\\] ]", ""));
                }
                Point point = new Point(i++, positionPoint, mark, -1, -1);
                pointList.add(point);
//                String key = values[mark].trim();
//                countMap.put(key, countMap.getOrDefault(key, 0) + 1);
                label_numbers.add(mark);
            }
            OutPut.k = label_numbers.size();
//            System.out.println("Value\tCount");
//            for (Map.Entry<String, Integer> entry : countMap.entrySet()) {
//                System.out.println(entry.getKey() + "\t" + entry.getValue());
//            }
        } catch (IOException e) {
            e.printStackTrace();
        }


        return pointList;
    }
}
