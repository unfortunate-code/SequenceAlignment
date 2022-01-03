// Author: Anirudh Alameluvari (alameluvarianirudh3@gmail.com)
// Date: 8th December 2021
import java.util.List;

public interface SequenceAlignment {

    /**
     * Aligns the given sequences with minimum cost.
     *
     * @param sequence1
     * @param sequence2
     */
    public void alignSequences(String sequence1, String sequence2);

    /**
     * Gets the final alignment.
     *
     * @return
     */
    public List<String> getAlignedSequences();

    /**
     * Gets the memory usage of the sequence alignment implementation in bytes.
     *
     * @return
     */
    public long getMemoryUsageInKBs();

    /**
     * Gets the time used by the sequence alignment implementation in nanoseconds.
     *
     * @return
     */
    public long getTimeInMillis();

    /**
     * Cost of alignment.
     *
     * @return
     */
    public int getAlignmentCost();
}