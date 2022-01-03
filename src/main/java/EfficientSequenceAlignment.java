// Author: Anirudh Alameluvari (alameluvarianirudh3@gmail.com)
// Date: 8th December 2021

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Does the space efficient sequence alignment.
 */
public class EfficientSequenceAlignment implements SequenceAlignment {
    private List<String> alignedSequences;
    private long memoryUsedInKBs;
    private long timeInMs;
    private int alignmentCost;
    private List<List<Integer>> path = new ArrayList<>();

    @Override
    public void alignSequences(String sequence1, String sequence2) {
        long start = System.currentTimeMillis();
        long memoryUsedAtStart = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
        for (int i = 0; i <= sequence1.length(); i++) {
            path.add(new ArrayList<>());
        }
        alignmentHelper(sequence1, sequence2);
        timeInMs = System.currentTimeMillis() - start;
        long memoryUsedAtEnd = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
        memoryUsedInKBs = (memoryUsedAtEnd - memoryUsedAtStart) / 1024;
    }

    /**
     * Manages the entire alignment logic and sets up the variables.
     *
     * @param X
     * @param Y
     */
    private void alignmentHelper(String X, String Y) {
        String Y1 = new StringBuilder(Y).reverse().toString();
        String X1 = new StringBuilder(X).reverse().toString();
        int[][] dpL = new int[Y.length() + 1][2];
        int[][] dpR = new int[Y1.length() + 1][2];
        divideAndConquerAlignment(X, Y, Y1, X1, dpL, dpR, 0, X.length() - 1, 0, Y.length() - 1);
        buildSequenceAndGetCost(path, X, Y);
    }

    /**
     * Main method for D&C alignment.
     *
     * @param X
     * @param Y
     * @param Y1
     * @param X1
     * @param dpL
     * @param dpR
     * @param lx
     * @param rx
     * @param ly
     * @param ry
     */
    private void divideAndConquerAlignment(String X, String Y, String Y1, String X1,
                                           int[][] dpL, int[][] dpR, int lx, int rx, int ly, int ry) {
        int m = rx - lx + 1;
        int n = ry - ly + 1;
        // Perform the basic alignment if one of the strings is less than 3 characters long.
        if (m <= 2 || n <= 2) {
            buildPathForBasic(basicAlignment(X, Y, lx, rx, ly, ry), X, Y, lx, rx, ly, ry);
            return;
        }
        // Get the cost array from forward and backward on left and right part of X respectively.
        spaceEfficientAlignment(X, Y, lx, (rx + lx) / 2, ly, ry, dpL);
        spaceEfficientAlignment(X1, Y1, X.length() - rx - 1, X.length() - (rx + lx) / 2 - 2,
                Y.length() - ry - 1, Y.length() - ly - 1, dpR);
        int q = 0;
        // Get the point in Y that minimizes the total cost.
        for (int i = 1; i <= n; i++) {
            if (dpL[q][1] + dpR[n - q][1] > dpL[i][1] + dpR[n - i][1]) {
                q = i;
            }
        }
        // Recursively call the same on the left and right parts.
        divideAndConquerAlignment(X, Y, Y1, X1, dpL, dpR, lx, (lx + rx) / 2, ly, ly + q - 1);
        divideAndConquerAlignment(X, Y, Y1, X1, dpL, dpR, (lx + rx) / 2 + 1, rx, ly + q, ry);
    }

    /**
     * Does the space efficient alignment.
     *
     * @param X
     * @param Y
     * @param lx
     * @param rx
     * @param ly
     * @param ry
     * @param dp
     */
    private void spaceEfficientAlignment(String X, String Y, int lx, int rx, int ly, int ry, int[][] dp) {
        int n = ry - ly + 1;
        int m = rx - lx + 1;
        for (int i = 0; i <= n; i++) {
            dp[i][0] = i * SequenceAlignmentParameters.GAP_PENALTY;
        }
        for (int j = 1; j <= m; j++) {
            dp[0][1] = j * SequenceAlignmentParameters.GAP_PENALTY;
            for (int i = 1; i <= n; i++) {
                dp[i][1] = dp[i - 1][0] +
                        SequenceAlignmentParameters.MISMATCH_COSTS[
                                SequenceAlignmentParameters.CHAR_MAP.get(Y.charAt(ly + i - 1))][
                                SequenceAlignmentParameters.CHAR_MAP.get(X.charAt(lx + j - 1))];
                dp[i][1] = Math.min(dp[i][1],
                        Math.min(dp[i - 1][1] + SequenceAlignmentParameters.GAP_PENALTY,
                        dp[i][0] + SequenceAlignmentParameters.GAP_PENALTY));
            }
            for (int i = 0; i <= n; i++) {
                dp[i][0] = dp[i][1];
            }
        }
    }

    /**
     * Does basic alignment and returns the dp array.
     *
     * @param X
     * @param Y
     * @param lx
     * @param rx
     * @param ly
     * @param ry
     * @return
     */
    private int[][] basicAlignment(String X, String Y, int lx, int rx, int ly, int ry) {
        int m = rx - lx + 1;
        int n = ry - ly + 1;
        int[][] dp = new int[n + 1][m + 1];
        for (int i  = 0; i <= n; i++) {
            dp[i][0] = i * SequenceAlignmentParameters.GAP_PENALTY;
        }
        for (int j = 0; j <= m; j++) {
            dp[0][j] = j * SequenceAlignmentParameters.GAP_PENALTY;
        }
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                dp[i][j] = dp[i - 1][j - 1] +
                        SequenceAlignmentParameters.MISMATCH_COSTS[
                                SequenceAlignmentParameters.CHAR_MAP.get(Y.charAt(ly + i - 1))][
                                SequenceAlignmentParameters.CHAR_MAP.get(X.charAt(lx + j - 1))];
                dp[i][j] = Math.min(Math.min(dp[i][j],
                        dp[i - 1][j] + SequenceAlignmentParameters.GAP_PENALTY),
                        dp[i][j - 1] + SequenceAlignmentParameters.GAP_PENALTY);
            }
        }
        return dp;
    }

    /**
     * Builds the path for the basic alignment.
     *
     * @param dp
     * @param X
     * @param Y
     * @param lx
     * @param rx
     * @param ly
     * @param ry
     */
    private void buildPathForBasic(int[][] dp, String X, String Y, int lx, int rx, int ly, int ry) {
        int i = ry - ly + 1;
        int j = rx - lx + 1;
        while (i >= 0 && j >= 0) {
            if (i == 0 && j == 0) break;
            int cost1 = Integer.MAX_VALUE;
            if (i > 0 && j > 0) {
                cost1 = dp[i - 1][j - 1] + SequenceAlignmentParameters.MISMATCH_COSTS[
                        SequenceAlignmentParameters.CHAR_MAP.get(Y.charAt(ly + i - 1))][
                        SequenceAlignmentParameters.CHAR_MAP.get(X.charAt(lx + j - 1))];
            }
            int cost2 = Integer.MAX_VALUE;
            if (i > 0) {
                cost2 = dp[i - 1][j] + SequenceAlignmentParameters.GAP_PENALTY;
            }
            int cost3 = Integer.MAX_VALUE;
            if (j > 0) {
                cost3 = dp[i][j - 1] + SequenceAlignmentParameters.GAP_PENALTY;
            }
            // This coordinate is passed by the dp array in building the efficient solution. Add it to the path.
            path.get(lx + j).add(ly + i);
            if (cost1 <= cost2 && cost1 <= cost3) {
                i--;
                j--;
            } else if (cost2 <= cost3) {
                i--;
            } else {
                j--;
            }
        }
    }

    /**
     * Builds the alignments and gets the alignment cost from the path.
     *
     * @param path
     * @param X
     * @param Y
     */
    private void buildSequenceAndGetCost(List<List<Integer>> path, String X, String Y) {
        StringBuilder alignment1 = new StringBuilder();
        StringBuilder alignment2 = new StringBuilder();
        List<Integer> start = path.get(0);
        Collections.sort(start);
        // 0th index of path corresponds to the entries in Y that need to be added before even the first index in X gets
        // added. Add _ for X and the corresponding index for Y.
        for (int y : start) {
            alignment1.append("_");
            alignment2.append(Y.charAt(y - 1));
            alignmentCost += SequenceAlignmentParameters.GAP_PENALTY;
        }
        int prevY = 0;
        if (start.size() > 0) {
            prevY = start.get(start.size() - 1);
        }
        for (int x = 1; x < path.size(); x++) {
            List<Integer> l = path.get(x);
            Collections.sort(l);
            int firstY = l.get(0);
            if (firstY != prevY) {
                alignment1.append(X.charAt(x - 1));
                alignment2.append(Y.charAt(firstY - 1));
                alignmentCost += SequenceAlignmentParameters.MISMATCH_COSTS[
                        SequenceAlignmentParameters.CHAR_MAP.get(Y.charAt(firstY - 1))][
                        SequenceAlignmentParameters.CHAR_MAP.get(X.charAt(x - 1))];
            } else {
                // If the next X gets the same Y as the last Y of previous X, this means no new Y was added to the path. Add _
                // to Y and the corresponding index of X here.
                alignment1.append(X.charAt(x - 1));
                alignment2.append("_");
                alignmentCost += SequenceAlignmentParameters.GAP_PENALTY;
            }
            for (int i = 1; i < l.size(); i++) {
                // If multiple values of Y are at the same X index, this means no new X is added for the subsequent Y
                // indices. Add _ for X and the corresponding index for Y.
                alignment1.append("_");
                alignment2.append(Y.charAt(l.get(i) - 1));
                alignmentCost += SequenceAlignmentParameters.GAP_PENALTY;
            }
            // Update prevY as the last Y of this X.
            prevY = l.get(l.size() - 1);
        }
        alignedSequences = new ArrayList<>();
        alignedSequences.add(alignment1.toString());
        alignedSequences.add(alignment2.toString());
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