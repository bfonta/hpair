package main

import (
	"fmt"
	"os"
	"io"
	"flag"
	"log"
	"errors"
	fp "path/filepath"
)

func check(e error) {
    if e != nil {
        panic(e)
    }
}

func main() {
	flag.Usage = func() {
	    fmt.Fprintf(os.Stderr, "Usage: go run create_files.go 1 2 3\n")
	}   
	flag.Parse()
	tail := flag.Args()
	
	base_dir := fp.Join(string(os.Getenv("HOME")), "hpair")
	out_dir  := fp.Join(base_dir, "inputs/")
	if _, err := os.Stat(out_dir); errors.Is(err, os.ErrNotExist) {
		os.Mkdir(out_dir, os.ModeDir)
	}

	// open template file
	template, err := os.Open("hpair.in")
    if err != nil {
        log.Fatal(err)
    }
    defer template.Close()
	
	for _, tri := range tail {
		fmt.Println(tri)
		dest, err := os.OpenFile(fp.Join(out_dir, "hpair_" + tri + ".txt"), os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0755)
		if err != nil {
			log.Fatal(err)
		}

		_, err = io.Copy(dest, template) // check first var for number of bytes copied
		check(err)

		err = dest.Sync()
		check(err)
		if err := dest.Close(); err != nil {
			log.Fatal(err)
		}

	}

	f, err := os.Create(fp.Join(out_dir, "hpair1.out"))
	check(err)
	defer f.Close()
}
