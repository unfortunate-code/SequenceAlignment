// Author: Anirudh Alameluvari (alameluvarianirudh3@gmail.com)
// Date: 8th December 2021
import java.util.ArrayList;
import java.util.List;

public class BasicSequenceAlignment implements SequenceAlignment {
    private List<String> alignedSequences;
    private long memoryUsedInKBs;
    private long timeInMs;
    private int alignmentCost;

    @Override
    public void alignSequences(String sequence1, String sequence2) {
        long start = System.currentTimeMillis();
        long memoryUsedAtStart = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
        alignmentHelper(sequence1, sequence2);
        timeInMs = System.currentTimeMillis() - start;
        long memoryUsedAtEnd = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
        memoryUsedInKBs = (memoryUsedAtEnd - memoryUsedAtStart) / 1024;
    }

    private void alignmentHelper(String sequence1, String sequence2) {
        int[][] dp;
        dp = new int[sequence1.length() + 1][sequence2.length() + 1];
        populateDPArray(sequence1, sequence2, dp);
        alignmentCost = dp[sequence1.length()][sequence2.length()];
        alignedSequences = new ArrayList<>();
        findAlignedSequences(sequence1, sequence2, dp);
    }

    /**
     * Backtracks to get the aligned sequences.
     *
     * @param sequence1
     * @param sequence2
     * @param dp
     */
    private void findAlignedSequences(String sequence1, String sequence2, int[][] dp) {
        StringBuilder aligned1 = new StringBuilder();
        StringBuilder aligned2 = new StringBuilder();
        int i = sequence1.length();
        int j = sequence2.length();
        while (i >= 0 && j >= 0) {
            if (i == 0 && j == 0) break;
            int cost1 = Integer.MAX_VALUE;
            // Cost for both Xi and Yj in the alignment
            if (i > 0 && j > 0) {
                 cost1 = dp[i - 1][j - 1] + SequenceAlignmentParameters.MISMATCH_COSTS[
                        SequenceAlignmentParameters.CHAR_MAP.get(sequence1.charAt(i - 1))][
                        SequenceAlignmentParameters.CHAR_MAP.get(sequence2.charAt(j - 1))];
            }
            int cost2 = Integer.MAX_VALUE;
            // Xi is not in the alignment.
            if (i > 0) {
                cost2 = dp[i - 1][j] + SequenceAlignmentParameters.GAP_PENALTY;
            }
            int cost3 = Integer.MAX_VALUE;
            // Yj is not in the alignment.
            if (j > 0) {
                cost3 = dp[i][j - 1] + SequenceAlignmentParameters.GAP_PENALTY;
            }
            // Find the minimum and create aligned sequences.
            if (cost1 <= cost2 && cost1 <= cost3) {
                aligned1.insert(0, sequence1.charAt(i - 1));
                aligned2.insert(0, sequence2.charAt(j - 1));
                i = i - 1;
                j = j - 1;
            } else if (cost2 <= cost3) {
                aligned1.insert(0, sequence1.charAt(i - 1));
                aligned2.insert(0, '_');
                i = i - 1;
            } else {
                aligned1.insert(0, '_');
                aligned2.insert(0, sequence2.charAt(j - 1));
                j = j - 1;
            }
        }
        alignedSequences.add(aligned1.toString());
        alignedSequences.add(aligned2.toString());
    }

    /**
     * Get costs bottom up.
     *
     * @param sequence1
     * @param sequence2
     * @param dp
     */
    private void populateDPArray(String sequence1, String sequence2, int[][] dp) {
        for (int i = 0; i <= sequence1.length(); i++) {
            dp[i][0] = i * SequenceAlignmentParameters.GAP_PENALTY;
        }
        for (int j = 0; j <= sequence2.length(); j++) {
            dp[0][j] = j * SequenceAlignmentParameters.GAP_PENALTY;
        }
        for (int i = 1; i <= sequence1.length(); i++) {
            for (int j = 1; j <= sequence2.length(); j++) {
                dp[i][j] = dp[i - 1][j - 1] +
                        SequenceAlignmentParameters.MISMATCH_COSTS[
                                SequenceAlignmentParameters.CHAR_MAP.get(sequence1.charAt(i - 1))][
                                SequenceAlignmentParameters.CHAR_MAP.get(sequence2.charAt(j - 1))];
                dp[i][j] = Math.min(Math.min(dp[i][j],
                        dp[i - 1][j] + SequenceAlignmentParameters.GAP_PENALTY),
                        dp[i][j - 1] + SequenceAlignmentParameters.GAP_PENALTY);
            }
        }
    }

    @Override
    public List<String> getAlignedSequences() {
        return alignedSequences;
    }

    @Override
    public long getMemoryUsageInKBs() {
        return memoryUsedInKBs;
    }

    @Override
    public long getTimeInMillis() {
        return timeInMs;
    }

    @Override
    public int getAlignmentCost() {
        return alignmentCost;
    }
}