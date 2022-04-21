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

	"github.com/spf13/cobra"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var ()

////////////////////////////////////////////////////////////////////////////////////////////////////

// accessionCmd represents the accession command
var accessionCmd = &cobra.Command{
	Use:   "accession",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("Accession called")
	},
}

func init() {
	rootCmd.AddCommand(accessionCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// accessionCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// accessionCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// package main

// ////////////////////////////////////////////////////////////////////////////////////////////////////

// import (
// 	"bufio"
// 	"log"
// 	"os"
// 	"strings"
// )

// ////////////////////////////////////////////////////////////////////////////////////////////////////

// // declarations
// var (
// 	readFile string = os.Args[1]
// 	fileOut  string = os.Args[2]
// )

// ////////////////////////////////////////////////////////////////////////////////////////////////////

// func main() {

// 	// execute logic
// 	proteinAccessionCollect(readFile)
// }

// ////////////////////////////////////////////////////////////////////////////////////////////////////

// // collect accessions
// func proteinAccessionCollect(readFile string) {

// 	// open an input file, exit on error
// 	inputFile, readErr := os.Open(readFile)
// 	if readErr != nil {
// 		log.Fatal("Error opening input file : ", readErr)
// 	}

// 	// scanner.Scan() advances to the next token returning false if an error was encountered
// 	scanner := bufio.NewScanner(inputFile)

// 	// iterate
// 	for scanner.Scan() {

// 		// collect pattern
// 		if strings.Contains(scanner.Text(), "protein_id") {

// 			records := strings.Split(scanner.Text(), "=")
// 			accession := strings.ReplaceAll(records[1], "\"", "")

// 			// write
// 			writeProtRecord(fileOut, accession)
// 		}
// 	}
// }

// ////////////////////////////////////////////////////////////////////////////////////////////////////

// // write records
// func writeProtRecord(fileOut, accession string) {

// 	// declare io
// 	f, err := os.OpenFile(fileOut, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0666)

// 	if err != nil {
// 		panic(err)
// 	}

// 	defer f.Close()

// 	// deslcare writer
// 	w := bufio.NewWriter(f)

// 	// writing
// 	_, err = w.WriteString(accession + "\n")

// 	if err != nil {
// 		panic(err)
// 	}

// 	// flush writer
// 	w.Flush()
// }

////////////////////////////////////////////////////////////////////////////////////////////////////
