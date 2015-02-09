***************
Getting Started
***************

VPN on Windows computers
========================
Several people have had problems using VPN while on Windows computers with
VirtualBox running LinuxMint OS.
Here are Erin Baschal's notes on Mike Campbell's solutions:

.. code-block:: bash

    How to fix internet/VPN issue in VirtualBox
    Open Network and Sharing Center
         Windows 8: right click on network icon at bottom of screen
    Change adapter settings (on left)
    Virtual Box Host Adapter (right click)
        Properties
        Uncheck TCP/IPv6
        Double click TCP/IPv4
        Use the following DNS server addresses:
            140.226.189.35
            132.194.70.65
        Advanced, DNS tab
            Append these DNS suffixes:
                ucdenver.pvt

