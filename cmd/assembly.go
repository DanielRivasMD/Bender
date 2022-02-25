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
	"math"
	"regexp"
	"strconv"

	"github.com/spf13/cobra"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var (
	assembly   string
	species    string
	outFile    string
	syncytin   identified // identified struct
	scaffoldID string
	startCoor  float64
	endCoor    float64
	hood       float64
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// syncytin features
type identified struct {
	scaffoldIdent string
	positionIdent position
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// positions
type position struct {
	startPos float64
	endPos   float64
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// annotations
type annotation struct {
	scaffoldIdent   string
	classIdent      string
	positionIdent   position
	scoreInt        int64
	strandSense     string
	attributeStruct attribute
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// attributes
type attribute struct {
	ID     string
	Name   string
	Alias  string
	Parent string
	Target string
	Note   string
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// TODO: documentation assembly
// assemblyCmd represents the assembly command
var assemblyCmd = &cobra.Command{
	Use:   "assembly",
	Short: "Handle assembly operations.",
	Long: `Handle assembly operations, such as:

description: parse assembly features.
loci: read diamond assembly search output & filter results.
search: perform similarity search.
sequence: extract sequences from assemblies.
synteny: .
`,

	Example: ``,
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	rootCmd.AddCommand(assemblyCmd)

	// persistent flags
	assemblyCmd.PersistentFlags().StringVarP(&assembly, "assembly", "a", "", "Assembly file")
	assemblyCmd.PersistentFlags().StringVarP(&species, "species", "p", "", "Species file")
	assemblyCmd.PersistentFlags().StringVarP(&outFile, "outfile", "o", "", "Out file. If empty it will be defined by input")

	assemblyCmd.PersistentFlags().StringVarP(&scaffoldID, "scaffold", "S", "", "Scaffold ID")
	assemblyCmd.PersistentFlags().Float64VarP(&startCoor, "start", "s", 0, "Start coordinate")
	assemblyCmd.PersistentFlags().Float64VarP(&endCoor, "end", "e", 0, "End coordinate")
	assemblyCmd.PersistentFlags().Float64VarP(&hood, "hood", "n", 0, "Neighborhood to look into")

}

////////////////////////////////////////////////////////////////////////////////////////////////////

// pass struct as reference to update
func (position *position) updateMinMax(num1, num2 float64) {
	position.startPos = math.Min(num1, num2)
	position.endPos = math.Max(num1, num2)
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// determine minimum / maximum
func (position *position) minMax(num1, num2 float64) {
	position.updateMinMax(num1, num2)
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// parse values & determine minimum / maximum
func (position *position) parseMinMax(str1, str2 string) {
	num1, _ := strconv.ParseFloat(str1, 64)
	num2, _ := strconv.ParseFloat(str2, 64)

	position.minMax(num1, num2)
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// define output file
func coordinateOut(readFile string) string {

	outFile = readFile
	reg := regexp.MustCompile(`HiC*`)
	res := reg.FindStringIndex(outFile)
	outFile = outFile[0:res[0]]

	// assembly directory
	outFile = outDir + "/" +
		outFile +
		syncytin.scaffoldIdent + "_" +
		strconv.FormatFloat(syncytin.positionIdent.startPos, 'f', 0, 64) + "_" +
		strconv.FormatFloat(syncytin.positionIdent.endPos, 'f', 0, 64)

	// return full path
	return outFile
}

////////////////////////////////////////////////////////////////////////////////////////////////////
