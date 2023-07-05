import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

public class Test {
    public static List<String> readSampleOutput(String fileName) throws IOException {
       return Files.lines(Path.of(fileName)).skip(1).map(s -> s.split("\t")[0]).collect(Collectors.toList());
    }

    public static void main(String[] args) throws IOException {
        List<String> outNodes = readSampleOutput("simul_exp_go_bp_ensembl_min50_max500.enrich.out");
    }
}
