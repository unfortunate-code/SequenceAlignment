// Author: Anirudh Alameluvari (alameluvarianirudh3@gmail.com)
// Date: 8th December 2021
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

/**
 * Entry point for sequence alignment.
 */
public class SequenceAlignmentMain {
    public static void main(String[] args) {
        List<String> sequences = parseInput();
        SequenceAlignment sequenceAlignment;
        if (args[0].equals(SequenceAlignmentParameters.RUN_BASIC_ALGORITHM)) {
            sequenceAlignment = new BasicSequenceAlignment();
        } else {
            sequenceAlignment = new EfficientSequenceAlignment();
        }
        sequenceAlignment.alignSequences(sequences.get(0), sequences.get(1));
        System.out.println(sequenceAlignment.getAlignedSequences().get(0));
        System.out.println(sequenceAlignment.getAlignedSequences().get(1));
        System.out.println(sequenceAlignment.getAlignmentCost());
        System.out.println(sequenceAlignment.getTimeInMillis());
        System.out.println(sequenceAlignment.getMemoryUsageInKBs());
    }

    private static List<String> parseInput() {
        Scanner scanner = new Scanner(System.in);
        String seed1 = scanner.nextLine();
        String s = scanner.nextLine();
        List<Integer> indices1 = new ArrayList<>();
        Integer num;
        while ((num = getNumber(s)) != null) {
            indices1.add(num);
            s = scanner.nextLine();
        }
        String seed2 = s;
        List<Integer> indices2 = new ArrayList<>();
        while (scanner.hasNextLine()) {
            indices2.add(getNumber(scanner.nextLine()));
        }
        List<String> sequences = new ArrayList<>();
        sequences.add(InputGenerator.generateInput(seed1, indices1));
        sequences.add(InputGenerator.generateInput(seed2, indices2));
        return sequences;
    }

    /**
     * Returns the number if the string is numeric. Otherwise, it returns null.
     *
     * @param s
     * @return
     */
    private static Integer getNumber(String s) {
        try {
            int num = Integer.parseInt(s);
            return num;
        } catch (NumberFormatException e) {
            return null;
        }
    }
}