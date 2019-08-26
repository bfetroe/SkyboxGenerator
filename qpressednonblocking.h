/* QPressedNonBlocking.h
 * Original author: Brandon Fetroe
 * Revised: 2018-09-10
 */


#ifndef QPRESSEDNONBLOCKING_H
#define QPRESSEDNONBLOCKING_H

#include <iostream>
#include <string>
#include <stdio.h>
#include <conio.h>
#include <sstream>

class qPressedNonBlocking
{
public:
    /* qPressed() returns true if q was pressed, false if it was not.
     * This method is meant to be used as a non-blocking way to detect if
     * the letter q was pressed (often while other work was being done).
     * Optionally, the whole string of keys that were pressed can be returned
     * by using qPressed(std::string& keysPressed).
     *
     * example usage if we want to return the string keysPressed:
     *
     * qPressedNonBlocking qPressedHandler;
     * std::string keysPressed;
     * for(int ii=0, ii<bigNumber, ii++){
     *     //...work being done...
     *     bool qWasPressed = qPressedHandler.qPressed(keysPressed);
     *     if(qWasPressed){
     *         cout << "Qutting program because 'q' was pressed"<<endl
     *              << "Full string entered that stopped this program = "
     *              << keysPressed <<endl;
     *         return 0;
     *     }
     *     //q wasn't pressed, so keep working...
     * }
     */

    bool qPressed(){
        std::string keys;
        return qPressed(keys);
    }

    bool qPressed(std::string& keysPressed){
        bool qWasPressed = false;
        keysPressed.clear();
        if(_kbhit()){
            std::stringstream ss;
            std::string currCharString;
            while(_kbhit()){
                ss.clear();
                currCharString.clear();
                char myChar = getch();
                ss << myChar;
                ss >> currCharString;
                keysPressed.append(currCharString);
                if(myChar == 'q'){
                    qWasPressed = true;
                }
            }
        }
        return qWasPressed;
    }
};

#endif // QPRESSEDNONBLOCKING_H
