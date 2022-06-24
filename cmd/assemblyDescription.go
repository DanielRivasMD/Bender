/*
Copyright © 2021 Daniel Rivas <danielrivasmd@gmail.com>

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
	"encoding/json"
	"io/ioutil"
	"log"
	"os"
	"strings"

	"github.com/spf13/cobra"
	"github.com/ttacon/chalk"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var ()

////////////////////////////////////////////////////////////////////////////////////////////////////

type Organism struct {
	Vernacular string
	Binominal  string
	FunFact    string
	Taxonomy   string
}

////////////////////////////////////////////////////////////////////////////////////////////////////

type ChromlengthAssembly struct {
	Karyotype         string
	ScaffoldN50       string
	NumberOfScaffolds string
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// data struct for the decoded data all fields must be exportable
type Data struct {
	Organism            Organism
	ChromlengthAssembly ChromlengthAssembly
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// descriptionCmd represents the description command
var descriptionCmd = &cobra.Command{
	Use:   "description",
	Short: "Collect assembly information.",
	Long: chalk.Green.Color("Daniel Rivas <danielrivasmd@gmail.com>") + `

Parse JSON file containing description of genome assemblies.

Collect information and concatenate.
`,

	Example: ``,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(κ *cobra.Command, args []string) {

		// execute logic
		describeAssembly(species)

	},
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	assemblyCmd.AddCommand(descriptionCmd)

	// flags

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func describeAssembly(readFile string) {

	// json file
	contentFile, ε := ioutil.ReadFile(inDir + "/" + readFile)
	if ε != nil {
		log.Fatal("Error when opening file: ", ε)
	}

	// unmarshall the data into features
	var features Data
	ε = json.Unmarshal(contentFile, &features)
	if ε != nil {
		log.Fatal("Error during Unmarshal(): ", ε)
	}

	// write extracted
	writeSpeciesFeatures(features)

}

////////////////////////////////////////////////////////////////////////////////////////////////////

// write positions
func writeSpeciesFeatures(features Data) {

	// declare io
	ƒ, ε := os.OpenFile(outDir+"/"+outFile, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0666)
	if ε != nil {
		panic(ε)
	}

	defer ƒ.Close()

	// declare writer
	ϖ := bufio.NewWriter(ƒ)

	// writing
	_, ε = ϖ.WriteString(
		features.Organism.Binominal + "," +
			strings.ReplaceAll(features.Organism.Vernacular, ",", "") + "," +
			features.ChromlengthAssembly.Karyotype + "," +
			strings.ReplaceAll(features.ChromlengthAssembly.ScaffoldN50, ",", "") + "," +
			strings.ReplaceAll(features.ChromlengthAssembly.NumberOfScaffolds, ",", "") +
			"\n")
	if ε != nil {
		panic(ε)
	}

	// flush writer
	ϖ.Flush()
}

////////////////////////////////////////////////////////////////////////////////////////////////////
