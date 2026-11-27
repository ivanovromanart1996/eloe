# EloE

EloE is a Java desktop application that calculates gene expression indices from annotated genomes in GenBank (`.gbk`) format. It provides both a graphical interface and a command-line mode.

EloE was written by V. S. Sokolov, B.S. Zuraev and S.A. Lashin. A. I. Klimenko added the command-line mode.

## Run EloE

1. Download `EloE-0.1.2.zip` from the **Releases** section and extract the entire archive.
2. Install Java 21 and check your installation with `java -version`. The graphical interface requires a desktop session.
3. Open a terminal in the extracted folder. Keep `EloE_desktop-0.1.2.jar` and the `lib` directory together.
4. Launch the graphical interface:

   ```sh
   java -jar EloE_desktop-0.1.2.jar
   ```

Create an input directory with one subdirectory per organism. Put one `.gbk` file in each subdirectory, for example `use/organism-1/genome.gbk`. Create a separate, empty output directory. In the application, select the parent input directory next to **Choose directory with organism/s**, select the output directory next to **Choose directory for results**, and click **Start**.

To run from a terminal, pass the paths to the input and output directories. Both directories must already exist. The following command works in PowerShell, macOS, and Linux:

```sh
java -cp EloE_desktop-0.1.2.jar org.eloe.console.ConsoleStarter "use" "res"
```

EloE writes `organisms_index.txt` and a subdirectory for each processed organism to the output directory. Choose a new, empty output directory for each run. EloE writes logs to the `log` directory in the current working directory.

## Build from source

Install JDK 21 and Maven 3.9. Run the following command from the directory containing `pom.xml`:

```sh
mvn clean package
```

Maven writes `target/EloE_desktop-0.1.2.jar` and copies the runtime dependencies to `target/lib/`. Package the JAR and the `lib` directory together when distributing the build. To launch the graphical interface from the project root, run `java -jar target/EloE_desktop-0.1.2.jar`.

## Known limitations

- The original documentation reported occasional hangs during stages 2 and 3 and a possible `NullPointerException` when calling `getOrganism`. These issues were reported for command-line mode and may also affect the graphical interface.
- Support for multiple GenBank records per organism was planned. The authors also reported cases of individual genes missing from the output.
- EloE uses the library versions listed in `pom.xml`. Third-party information is included in the dependency JARs where provided.
