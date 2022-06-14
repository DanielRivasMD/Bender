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
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
const readmeя = `README`
const assemblyя = `HiC.fasta.gz`
const annotя = `fasta_v2.functional.gff3.gz`

var (
	// declare regex
	readmeρ   = regexp.MustCompile(readmeя)
	assemblyρ = regexp.MustCompile(assemblyя)
	annotρ    = regexp.MustCompile(annotя)

	// declare switches
	assemblyϟ bool
	annotϟ    bool

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
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(κ *cobra.Command, args []string) {

		// execute logic
		assemblyFilter(inDir + "/" + species)

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
	inputFile, ε := os.Open(readFile)
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
		case readmeρ.MatchString(records[0]):
			readmeLink = records[1]

		case assemblyρ.MatchString(records[0]):
			assemblyϟ = true
			assemblyID = records[0]
			assemblyLink = records[1]

		case annotρ.MatchString(records[0]):
			annotϟ = true
			annotID = records[0]
			annotLink = records[1]
		}

	}

	if assemblyϟ && annotϟ {
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
