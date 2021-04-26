////////////////////////////////////////////////////////////////////////////////////////////////////

package cmd

////////////////////////////////////////////////////////////////////////////////////////////////////

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

type ratio struct {
	numerator   float64
	denominator float64
}

type blast struct {
	identities ratio
	positives  ratio
	gaps       ratio
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func blastFilter(preffixFile string) {

	// open an input file, exit on error
	inputFile, readErr := os.Open(preffixFile + ".out")
	if readErr != nil {
		log.Fatal("Error opening input file : ", readErr)
	}

	//scanner.Scan() advances to the next token returning false if an error was encountered
	scanner := bufio.NewScanner(inputFile)
	for scanner.Scan() {
		records := strings.Split(scanner.Text(), " ")

		// TODO: split scanner from filter

		// filter results
		// TODO: record location & identity
		if len(records) > 1 {
			var blast blast
			if strings.Contains(records[1], "Identities") {

				blast.identities.numerator, _ = strconv.ParseFloat(strings.Split(records[3], "/")[0], 64)
				blast.identities.denominator, _ = strconv.ParseFloat(strings.Split(records[3], "/")[1], 64)

				blast.positives.numerator, _ = strconv.ParseFloat(strings.Split(records[7], "/")[0], 64)
				blast.positives.denominator, _ = strconv.ParseFloat(strings.Split(records[7], "/")[1], 64)

				blast.gaps.numerator, _ = strconv.ParseFloat(strings.Split(records[11], "/")[0], 64)
				blast.gaps.denominator, _ = strconv.ParseFloat(strings.Split(records[11], "/")[1], 64)

				// if blast.identities.denominator > 400 && (blast.identities.numerator/blast.identities.denominator) > 0.8 {
				if blast.identities.denominator > 400 && (blast.positives.numerator/blast.positives.denominator) > 0.8 {

					fmt.Println(blast)
					text := fmt.Sprintf("%g | %g\n", blast.identities.numerator, blast.identities.denominator)

					// TODO: write to file
					fileOut := preffixFile + ".filtered.out"
					f, err := os.OpenFile(fileOut, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0666)

					if err != nil {
						panic(err)
					}

					defer f.Close()

					w := bufio.NewWriter(f)
					_, err = w.WriteString(text)
					if err != nil {
						panic(err)
					}

					w.Flush()

				}

			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////
