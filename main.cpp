/*
sortContigs uses a binary tree to input DNA sequences from a fasta file
it then outputs them to the designated output file, sorted by sequence length.

Copyright 2010 Mark Sechter
*/

#include <iostream>
#include <fstream>
#include <pcrecpp.h>
#include <stdlib.h>

using namespace std;

struct node
{
    int key_value;
    string record_value;
    node *left;
    node *right;
};

class btree
{
    public:
      btree();
      ~btree();

      void insert(int key, string record);
      node *search(int key);
      void destroy_tree();
      void outputTree(ofstream & writeFile);

    private:
      void destroy_tree(node *leaf);
      void insert(int key, string record, node *leaf);
      node *search(int key, node *leaf);
      void outputTree(node *leaf, ofstream & writeFile);

      node *root;
      string contigRecord;

};

btree::btree()
{
    root=NULL;
}

void btree::destroy_tree(node *leaf)
{
    if(leaf!=NULL)
    {
        destroy_tree(leaf->left);
        destroy_tree(leaf->right);
        delete leaf;
    }
}

void btree::insert(int key, string record, node *leaf)
{
    if(key < leaf->key_value)
    {
        if(leaf->left!=NULL)
          insert(key, record, leaf->left);
        else
        {
            leaf->left=new node;
            leaf->left->key_value=key;
            leaf->left->record_value=record;
            leaf->left->left=NULL;
            leaf->left->right=NULL;
        }
    }
    else if(key >= leaf->key_value)
    {
        if(leaf->right!=NULL)
          insert(key, record, leaf->right);
        else
        {
            leaf->right=new node;
            leaf->right->key_value=key;
            leaf->right->record_value=record;
            leaf->right->left=NULL;
            leaf->right->right=NULL;
        }
    }
}

void btree::insert(int key, string record)
{
    if(root!=NULL)
      insert(key, record, root);
    else
    {
        root=new node;
        root->key_value=key;
        root->record_value=record;
        root->left=NULL;
        root->right=NULL;
    }
}

btree::~btree()
{
    destroy_tree();
}

node *btree::search(int key, node *leaf)
{
    if(leaf!=NULL)
    {
        if(key==leaf->key_value)
          return leaf;
        if(key < leaf->key_value)
          return search(key, leaf->left);
        else
          return search(key, leaf->right);
    }
    else return NULL;
}

node *btree::search(int key)
{
    return search(key, root);
}

void btree::destroy_tree()
{
    destroy_tree(root);
}

void btree::outputTree(ofstream & writeFile)  //public function to output tree
{
    if(root!=NULL) //if the root node has a value
    {
      outputTree(root, writeFile);  //call the private outputTree function on it
    }
}

void btree::outputTree(node *leaf, ofstream & writeFile)
{
    if(leaf->left!=NULL)  //if the node has a left child, try to output its contig first
      outputTree(leaf->left, writeFile);

    writeFile << leaf->record_value << endl; //then output the current node's contig to the designated outfile

    if(leaf->right!=NULL) //then, if the node has a right child, try to output that node's contig
      outputTree(leaf->right, writeFile);
}

int main(int argc, char *argv[])  //accept command line arguments
{
    btree contigTree;  //declare a tre that will hold the unordered contigs
    char getch;  //a char to hold chars while reading the file
    string line;  //a string to hold each 2-line read as it is assembled with chars from getch
    string keyvalStr;  //a string to hold the contig length as it is extracted with regex
    int keyval;  //an int to hold the contig length extracted from the first line of each read
    bool isTwoLines = false;  //a bool to ignore 1 \n char (every read is 2 lines long)

    //should i use const_cast ??
    const char* fastaIn = argv[1];  // a pointer to a char array where we have the input fasta file name (cmd line arg 1)
    const char* fastaOut = argv[2]; // a pointer to a char array wher we have the output fasta file name (cmd line arg 2)

    if (argc != 3)
    {
      cout << "sortContigs is a utility to sort fasta files by contig length\n";
      cout << "Usage: " << argv[0] << " inputFileName outputFileName\n";
      cout << "Warning: Do not designate the same input and outpfile as your data will be lost.";
      return 0;
    }

    if (strcmp(argv[1], argv[2]) == 0)
    {
        cout << "inputFileName cannot be the same as outputFileName.\n";
        return 0;
    }


    ofstream outFile (fastaOut); //create an ofstream using the designated output file name
    if (!outFile.is_open())
    {
        cout << "Cannot write to file: " << fastaOut << endl;
        return 0;
    }

    ifstream inFile ( fastaIn );  //an input file stream for the input file
    if ( !inFile.is_open() )  //check if the file was opened. if so, continue
    {
      cout << "Could not open file: " << fastaIn << endl; // if not, notify user.
      return 0;  // and exit
    }
    else {
        cout << "Loading binary tree with contigs from file " << fastaIn << ".\n";
        while (!inFile.eof())  //check if the file has more data.  if so, continue.
        {
            inFile.get(getch);  //read the next char in the file into 'getch'
            if (getch=='\r') continue;  //ignore carriage return characters
            if (getch=='\n' || inFile.eof())  //if the char is a newline char, continue.
            {
                if (isTwoLines)  //if the line already has a newline char, continue.
                {
                    keyvalStr = line;  //copy the line into keyvalStr
                    pcrecpp::RE(">\\d+\\s(\\d+).+\\n.+").Replace("\\1", &keyvalStr);  //cut out all but the contig length from keyvalStr
                    keyval = atoi(keyvalStr.c_str());  //convert the contig length to an integer
                    contigTree.insert(keyval, line); //insert the contig length and it's associated line into the tree!

                    line=""; //clear the line
                    isTwoLines=false;  //indicate that line is now only 1 line long again.
                }
                else  // if the line is only 1 line long, continue
                {
                    line+=getch;  //add the newline char to the current line
                    isTwoLines=true;  //set the line to now be 2 lines long
                }
            }
            else  //if the char is not a newline char, continue
            {
                line+=getch; //add the non-newline char on to the current line
            }
        }

        cout << "Done.\nOutputting contigs from binary tree to file " << fastaOut << ".\n";
        contigTree.outputTree(outFile); //output all of the contigs now sorted by their length!
        cout << "Done.\n";

        return 0;
    }
}
