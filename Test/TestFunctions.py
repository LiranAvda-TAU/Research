from random import *


class TestFunctions:

    def __init__(self, function: str, dictionary: dict = {}, lst: list = [], result: bool = 0):
        self.function = function
        self.dictionary = dictionary
        self.lst = lst
        self.result = result
        print("-------------Test for " + self.function + " function has started-------------")

    def testFinisher(self):
        print("-------------Test for " + self.function + " function has finished-------------")
        print("\n")

    def checkBooleanResult(self):
        if self.result:
            print("Hooray! The function: " + self.function + " works!")
        else:
            print("Something is wrong with " + self.function + ". Try again :(")

    def checkSize(self, incorrectSize=0):
        if len(self.dictionary) != incorrectSize or len(self.lst) != incorrectSize:
            print("Number of items is: " + str(max(len(self.dictionary), len(self.lst))))
            self.result = True
        else:
            print("no luck")
            self.result = False

    def printFirstLinesInDict(self, x: int):
        size = len(self.dictionary)
        if size < x:
            x = size
            print("dictionary has only " + str(size) + " items")
        lst = list(self.dictionary.items())
        print(str(x) + " first lines are: ")
        for i in range(x):
            print(lst[i])

    def printFirstLinesInLst(self, x: int):
        size = len(self.lst)
        if size < x:
            x = size
            print("list has only " + str(size) + " items")
        print(str(x) + " first lines are: ")
        for i in range(x):
            print(self.lst[i])

    def printRandomLinesInDict(self, x: int):
        size = len(self.dictionary)
        if size < x:
            x = size
            print("dictionary has only " + str(size) + " items")
        items = list(self.dictionary.items())
        print(str(x) + " random lines are: ")
        samples = sample(items, x)  # Pick x random items from the list
        for s in samples:
            print(s)

    def printRandomLinesInLst(self, x: int):
        size = len(self.lst)
        if size < x:
            x = size
            print("list has only " + str(size) + " items")
        print(str(x) + " random lines are: ")
        samples = sample(self.lst, x)  # Pick x random items from the list
        for s in samples:
            print(s)


