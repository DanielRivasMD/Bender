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
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"strings"

	"github.com/atrox/homedir"
	"github.com/labstack/gommon/color"
	"github.com/spf13/cobra"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

var (
	assembly    string
	assemblyDir string
	species     string
	library     string
	libraryDir  string
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// AssemblySearchCmd represents the AssemblySearch command
var AssemblySearchCmd = &cobra.Command{
	Use:   "AssemblySearch",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	ValidArgs: []string{"blast", "diamond", "last"},
	Args:      cobra.ExactValidArgs(1),

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// Find home directory.
		home, errHomedir := homedir.Dir()
		if errHomedir != nil {
			fmt.Println(errHomedir)
			os.Exit(1)
		}

		// trim suffixes
		library = strings.TrimSuffix(library, ".fasta")
		assembly = strings.TrimSuffix(assembly, ".fasta.gz")

		// lineBreaks
		lineBreaks()

		// buffers
		var stdout bytes.Buffer
		var stderr bytes.Buffer

		searchType := ""
		switch args[0] {
		case "blast":
			searchType = "GenomeBlast.sh"
		case "diamond":
			searchType = "GenomeDiamond.sh"
		case "last":
			searchType = "GenomeLast.sh"
		}

		// shell call
		commd := home + "/bin/goTools/sh/" + searchType
		shCmd := exec.Command(commd, species, assembly, assemblyDir, library, libraryDir, outDir)

		// run
		shCmd.Stdout = &stdout
		shCmd.Stderr = &stderr
		_ = shCmd.Run()

		// stdout
		color.Println(color.Cyan(stdout.String(), color.B))

		// stderr
		if stderr.String() != "" {
			color.Println(color.Red(stderr.String(), color.B))
		}

		// lineBreaks
		lineBreaks()

	},
	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	rootCmd.AddCommand(AssemblySearchCmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags
	AssemblySearchCmd.Flags().StringVarP(&species, "species", "s", "", "Species")
	AssemblySearchCmd.Flags().StringVarP(&assembly, "assembly", "a", "", "Assembly to search")
	AssemblySearchCmd.Flags().StringVarP(&assemblyDir, "assemblyDir", "A", "", "Assembly directory")
	AssemblySearchCmd.Flags().StringVarP(&library, "library", "l", "", "Library to search against")
	AssemblySearchCmd.Flags().StringVarP(&libraryDir, "libraryDir", "L", "", "Library directory")

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////
