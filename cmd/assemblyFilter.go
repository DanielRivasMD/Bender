/*
Copyright © 2022 Daniel Rivas <danielrivasmd@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/
package cmd

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"regexp"
	"strings"

	"github.com/spf13/cobra"
	"github.com/ttacon/chalk"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
const rmя = `README`
const asя = `HiC.fasta.gz`
const anя = `fasta_v2.functional.gff3.gz`

var (
	// declare regex
	rmρ = regexp.MustCompile(rmя)
	asρ = regexp.MustCompile(asя)
	anρ = regexp.MustCompile(anя)

	// declare switches
	asϟ bool
	anϟ bool

	// declare placeholders
	readmeLink   string
	assemblyID   string
	annotID      string
	assemblyLink string
	annotLink    string
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// filterCmd represents the filter command
var filterCmd = &cobra.Command{
	Use:   "filter",
	Short: "Filter assembly list.",
	Long: chalk.Green.Color("Daniel Rivas <danielrivasmd@gmail.com>") + `

Filter assembly list based on the presence of chromosome length assembly & annotation.
`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(κ *cobra.Command, args []string) {

		// execute logic
		assemblyFilter(species)

	},
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	assemblyCmd.AddCommand(filterCmd)

	// flags

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func assemblyFilter(readFile string) {

	// open an input file, exit on error
	inputFile, ε := os.Open(inDir + "/" + readFile)
	if ε != nil {
		log.Fatal("Error opening input file : ", ε)
	}

	// scanner.Scan() advances to the next token returning false if an error was encountered
	scanner := bufio.NewScanner(inputFile)

	for scanner.Scan() {

		// comma separated records
		records := strings.Split(scanner.Text(), ",")

		// match records
		switch {
		case rmρ.MatchString(records[0]):
			readmeLink = records[1]

		case asρ.MatchString(records[0]):
			asϟ = true
			assemblyID = records[0]
			assemblyLink = records[1]

		case anρ.MatchString(records[0]):
			anϟ = true
			annotID = records[0]
			annotLink = records[1]
		}

	}

	// if asϟ && anϟ {
	if asϟ {
		fmt.Println(
			strings.Replace(species, ".csv", "", -1) + "," +
				assemblyID + "," +
				annotID + "," +
				readmeLink + "," +
				assemblyLink + "," +
				annotLink,
		)
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////
