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
const readme = `README`
const assemblyHiC = `HiC.fasta.gz`
const annotationGFF3 = `fasta_v2.functional.gff3.gz`

var (
	// declare regex
	я_readme = regexp.MustCompile(readme)
	я_assembly = regexp.MustCompile(assemblyHiC)
	я_annotation = regexp.MustCompile(annotationGFF3)

	// declare switches
	ϙ_assembly bool
	ϙ_annotation bool

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

// filter assemblies
func assemblyFilter(ƒ string) {

	// open an input file, exit on error
	inputFile, ε := os.Open(inDir + "/" + ƒ)
	if ε != nil {
		log.Fatal("Error opening input file : ", ε)
	}

	// scanner.Scan() advances to the next token returning false if an error was encountered
	scanner := bufio.NewScanner(inputFile)

	// iterate by lines
	for scanner.Scan() {

		// comma separated records
		records := strings.Split(scanner.Text(), ",")

		// match records
		switch {
		case я_readme.MatchString(records[0]):
			readmeLink = records[1]

		case я_assembly.MatchString(records[0]):
			ϙ_assembly = true
			assemblyID = records[0]
			assemblyLink = records[1]

		case я_annotation.MatchString(records[0]):
			ϙ_annotation = true
			annotID = records[0]
			annotLink = records[1]
		}

	}

	// if asϟ && anϟ {
	if ϙ_assembly {
		writeAssemblies()
	} else {
		fmt.Println(ƒ)
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////

// write assemblies
func writeAssemblies() {

	// declare io
	ƒ, ε := os.OpenFile(outDir+"/"+outFile, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0666)
	if ε != nil {
		panic(ε)
	}

	defer ƒ.Close()

	// declare writer
	writer := bufio.NewWriter(ƒ)

	// writing
	_, ε = writer.WriteString(

		strings.Replace(species, ".csv", "", -1) + "," +
			assemblyID + "," +
			annotID + "," +
			readmeLink + "," +
			assemblyLink + "," +
			annotLink +
			"\n")
	if ε != nil {
		panic(ε)
	}

	// flush writer
	writer.Flush()

}

////////////////////////////////////////////////////////////////////////////////////////////////////
