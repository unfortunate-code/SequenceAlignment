// Author: Anirudh Alameluvari (alameluvarianirudh3@gmail.com)
// Date: 8th December 2021
import java.util.HashMap;
import java.util.Map;

/**
 * Holds the hardcoded constants for sequence alignment.
 */
public class SequenceAlignmentParameters {

    public static final int GAP_PENALTY = 30;
    public static final int[][] MISMATCH_COSTS = {{0, 110, 48, 94}, {110, 0, 118, 48},
            {48, 118, 0, 110}, {94, 48, 110, 0}};
    public static final Map<Character, Integer> CHAR_MAP = new HashMap<>() {{
        put('A', 0);
        put('C', 1);
        put('G', 2);
        put('T', 3);
    }};
    public static final String RUN_BASIC_ALGORITHM = "basic";
}