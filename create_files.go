package main

import (
	"fmt"
	"os"
	"flag"
	"errors"
	"strings"
	io "io/ioutil"
	fp "path/filepath"
)

func check(e error) {
    if e != nil {
        panic(e)
		// log.Fatalln(err)
    }
}

func create_dir(d string) {
	if _, err := os.Stat(d); errors.Is(err, os.ErrNotExist) {
		os.Mkdir(d, 0700) //os.Mkdir(d, os.ModeDir)
	}
}

func main() {
	flag.Usage = func() {
	    fmt.Fprintf(os.Stderr, "Usage: go run create_files.go 1 2 3\n")
	}   
	flag.Parse()
	tail := flag.Args()

	base_dir := fp.Join(string(os.Getenv("HOME")), "hpair")
	out_dir := fp.Join(base_dir, "inputs/") // the inputs folder represents the output of this script
	create_dir(out_dir)
	create_dir(fp.Join(base_dir, "outputs/"))

	// open template file
	template, err := io.ReadFile("hpair.template.in")
	check(err)
    // defer template.Close()

	// create input files
	for _, tri := range tail {
		var tri2 string;
		if strings.Contains(tri, "m") {
			tri2 = strings.Replace(tri, "m", "-", -1)
		} else {
			tri2 = tri
		}
		fmt.Println(tri + " " + tri2)
		
		template_lines := strings.Split(string(template), "\n")
        for i, line := range template_lines {
			if strings.Contains(line, "C_3") {
				template_lines[i] = "C_3      =  " + tri2 + ".0D0";
			}
        }

		output_lines := strings.Join(template_lines, "\n")
        err = io.WriteFile(fp.Join(out_dir, "hpair_" + tri + ".in"), []byte(output_lines), 0644)
		check(err)
	}

	
}
