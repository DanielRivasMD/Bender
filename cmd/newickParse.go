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
	"fmt"
	"os"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"

	"github.com/spf13/cobra"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

var (
	newickFile string
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// newickParseCmd represents the newickParse command
var newickParseCmd = &cobra.Command{
	Use:   "newickParse",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// execute logic
		newickParse(inDir + "/" + newickFile)

	},

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	rootCmd.AddCommand(newickParseCmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags
	newickParseCmd.Flags().StringVarP(&newickFile, "newick", "n", "", "Newick file to parse")

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

// parse newick file
func newickParse(readFile string) {
	fmt.Println(readFile)

	var t *tree.Tree
	var err error
	var f *os.File
	if f, err = os.Open(readFile); err != nil {
		panic(err)
	}
	t, err = newick.NewParser(f).Parse()
	if err != nil {
		panic(err)
	}
	fmt.Println(t.Newick())

}

////////////////////////////////////////////////////////////////////////////////////////////////////
