// Author: Anirudh Alameluvari (alameluvarianirudh3@gmail.com)
// Date: 8th December 2021
import java.util.List;

public class InputGenerator {

    /**
     * Generates the input sequence from the user input.
     *
     * @param seed Seed given by the user.
     * @param indices The indices at which to duplicate.
     * @return
     */
    public static String generateInput(String seed, List<Integer> indices) {
        StringBuilder inputBuilder = new StringBuilder();
        inputBuilder.append(seed);
        for (int index : indices) {
            String s = inputBuilder.toString();
            inputBuilder.insert(index + 1, s);
        }
        return inputBuilder.toString();
    }
}