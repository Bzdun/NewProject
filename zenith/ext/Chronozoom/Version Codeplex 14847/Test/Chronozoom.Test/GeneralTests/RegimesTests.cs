﻿using System;
using System.Text;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using OpenQA.Selenium;
using OpenQA.Selenium.Interactions;
using OpenQA.Selenium.Interactions.Internal;
using Chronozoom.Test.Auxiliary;
using Chronozoom.Test.Components;
using System.Drawing.Imaging;
using System.Drawing;
using System.Threading;
using Chronozoom.Test.JsTypes;

namespace Chronozoom.Test.GeneralTests
{
    [TestClass]
    [TestPage(CzCommon.CzBetaStartPage)]
    public abstract class RegimesTests : CzTestBase
    {
        private VirtualCanvasComponent vcPageObj;
        private ActionsExtension action;

        [TestInitialize]
        public void TestInitialize()
        {
            GoToUrl();
            vcPageObj = new VirtualCanvasComponent(Driver);
            vcPageObj.WaitContentLoading();
        }

        [TestCleanup]
        public void TestCleanup()
        {
            if (action != null)
            {
                action.SetDefault();
            }
        }

        [TestMethod]
        public void TestRegimes_ClickHumanityLink_HumanityTimelineIsVisible()
        {
            Point linkPosition = vcPageObj.HumanityLink.Location;
            JsVisible visibleBefore = vcPageObj.GetViewport();

            action = new ActionsExtension(Driver);

            // Click on the regime link.
            action.MoveByOffset(linkPosition.X + 5, linkPosition.Y + 2).Click().Perform();

            vcPageObj.WaitAnimation();

            JsVisible visibleAfter = vcPageObj.GetViewport();
            JsTimeline timeline = GetHumanityTimeline();
            double width = timeline.right - timeline.left;
            JsCoordinates timelineCenter = new JsCoordinates(timeline.left + width / 2, timeline.y + timeline.height / 2);

            Assert.IsTrue(IsTimlineInFullSize(visibleAfter, timeline));
            Assert.AreEqual(visibleAfter.CenterX, timelineCenter.X, 1);
            Assert.AreEqual(visibleAfter.CenterY, timelineCenter.Y, 1);
            Assert.AreNotEqual(visibleBefore, visibleAfter);
        }

        [TestMethod]
        public void TestRegimes_ClickPrehistoryLink_PrehistoryTimelineIsVisible()
        {
            Point linkPosition = vcPageObj.PrehistoryLink.Location;
            JsVisible visibleBefore = vcPageObj.GetViewport();

            action = new ActionsExtension(Driver);

            // Click on the regime link.
            action.MoveByOffset(linkPosition.X + 5, linkPosition.Y + 2).Click().Perform();

            vcPageObj.WaitAnimation();

            JsVisible visibleAfter = vcPageObj.GetViewport();
            JsTimeline timeline = GetPrehistoryTimeline();
            double width = timeline.right - timeline.left;
            JsCoordinates timelineCenter = new JsCoordinates(timeline.left + width / 2, timeline.y + timeline.height / 2);

            Assert.IsTrue(IsTimlineInFullSize(visibleAfter, timeline));
            Assert.AreEqual(visibleAfter.CenterX, timelineCenter.X, 1);
            Assert.AreEqual(visibleAfter.CenterY, timelineCenter.Y, 1);
            Assert.AreNotEqual(visibleBefore, visibleAfter);
        }

        [TestMethod]
        public void TestRegimes_ClickLifeLink_LifeTimelineIsVisible()
        {
            Point linkPosition = vcPageObj.LifeLink.Location;
            JsVisible visibleBefore = vcPageObj.GetViewport();

            action = new ActionsExtension(Driver);

            // Click on the regime link.
            action.MoveByOffset(linkPosition.X + 5, linkPosition.Y + 2).Click().Perform();

            vcPageObj.WaitAnimation();

            JsVisible visibleAfter = vcPageObj.GetViewport();
            JsTimeline timeline = GetLifeTimeline();
            double width = timeline.right - timeline.left;
            JsCoordinates timelineCenter = new JsCoordinates(timeline.left + width / 2, timeline.y + timeline.height / 2);

            Assert.IsTrue(IsTimlineInFullSize(visibleAfter, timeline));
            Assert.AreEqual(visibleAfter.CenterX, timelineCenter.X, 1);
            Assert.AreEqual(visibleAfter.CenterY, timelineCenter.Y, 1);
            Assert.AreNotEqual(visibleBefore, visibleAfter);
        }

        [TestMethod]
        public void TestRegimes_ClickEarthLink_EarthTimelineIsVisible()
        {
            Point linkPosition = vcPageObj.EarthLink.Location;
            JsVisible visibleBefore = vcPageObj.GetViewport();

            action = new ActionsExtension(Driver);

            // Click on the regime link.
            action.MoveByOffset(linkPosition.X + 5, linkPosition.Y + 2).Click().Perform();

            vcPageObj.WaitAnimation();

            JsVisible visibleAfter = vcPageObj.GetViewport();
            JsTimeline timeline = GetEarthTimeline();
            double width = timeline.right - timeline.left;
            JsCoordinates timelineCenter = new JsCoordinates(timeline.left + width / 2, timeline.y + timeline.height / 2);

            Assert.IsTrue(IsTimlineInFullSize(visibleAfter, timeline));
            Assert.AreEqual(visibleAfter.CenterX, timelineCenter.X, 1);
            Assert.AreEqual(visibleAfter.CenterY, timelineCenter.Y, 1);
            Assert.AreNotEqual(visibleBefore, visibleAfter);
        }

        [TestMethod]
        public void TestRegimes_ClickCosmosLink_CosmosTimelineIsVisible()
        {
            Point cosmosLinkPosition = vcPageObj.CosmosLink.Location;
            Point earthLinkPosition = vcPageObj.EarthLink.Location;

            action = new ActionsExtension(Driver);

            // Click on the earth regime link.
            action.MoveByOffset(earthLinkPosition.X + 5, earthLinkPosition.Y + 2).Click().Perform();

            vcPageObj.WaitAnimation();

            action.SetDefault();

            JsVisible visibleBefore = vcPageObj.GetViewport();

            // Click on the cosmos regime link.
            action.MoveByOffset(cosmosLinkPosition.X + 5, cosmosLinkPosition.Y + 2).Click().Perform();

            vcPageObj.WaitAnimation();

            JsVisible visibleAfter = vcPageObj.GetViewport();
            JsTimeline timeline = GetCosmosTimeline();
            double width = timeline.right - timeline.left;
            JsCoordinates timelineCenter = new JsCoordinates(timeline.left + width / 2, timeline.y + timeline.height / 2);

            Assert.IsTrue(IsTimlineInFullSize(visibleAfter, timeline));
            Assert.AreEqual(visibleAfter.CenterX, timelineCenter.X, 1);
            Assert.AreEqual(visibleAfter.CenterY, timelineCenter.Y, 1);
            Assert.AreNotEqual(visibleBefore, visibleAfter);
        }

        public bool IsTimlineInFullSize(JsVisible visible, JsTimeline timeline)
        {
            Size vcSize = vcPageObj.VirtualCanvas.Size;
            double vwidth = timeline.right - timeline.left;
            double vheight = timeline.height;
            double pwidth = vwidth / visible.Scale;
            double pheight = vheight / visible.Scale;
            double timelineAspectRatio = pwidth / pheight;
            double canvasAspectRatio = (double)vcSize.Width / vcSize.Height;
            double ratio = 0;

            if (timelineAspectRatio >= canvasAspectRatio)
            {
                ratio = pwidth / vcSize.Width;
            }
            else
            {
                ratio = pheight / vcSize.Height;
            }

            // NOTE: It's possible to regulate this condition. The term of "full size" is documented.
            return ratio > 0.9 && ratio <= 1;
        }

        // NOTE: These methods get only important properties of timeline.
        public JsTimeline GetCosmosTimeline()
        {
            JsTimeline tl = new JsTimeline();

            tl.left = ExecuteScriptGetNumber("return content.d[0].left;");
            tl.right = ExecuteScriptGetNumber("return content.d[0].right;");
            tl.height = ExecuteScriptGetNumber("return content.d[0].height;");
            tl.realHeight = ExecuteScriptGetNumber("return content.d[0].realHeight;");
            tl.y = ExecuteScriptGetNumber("return content.d[0].y;");
            tl.realY = ExecuteScriptGetNumber("return content.d[0].realY;");

            return tl;
        }

        public JsTimeline GetEarthTimeline()
        {
            JsTimeline tl = new JsTimeline();

            tl.left = ExecuteScriptGetNumber("return FindChildTimeline(content.d[0], earthTimelineID).left;");
            tl.right = ExecuteScriptGetNumber("return FindChildTimeline(content.d[0], earthTimelineID).right;");
            tl.height = ExecuteScriptGetNumber("return FindChildTimeline(content.d[0], earthTimelineID).height;");
            tl.realHeight = ExecuteScriptGetNumber("return FindChildTimeline(content.d[0], earthTimelineID).realHeight;");
            tl.y = ExecuteScriptGetNumber("return FindChildTimeline(content.d[0], earthTimelineID).y;");
            tl.realY = ExecuteScriptGetNumber("return FindChildTimeline(content.d[0], earthTimelineID).realY;");

            return tl;
        }

        public JsTimeline GetLifeTimeline()
        {
            JsTimeline tl = new JsTimeline();
            string script = "return FindChildTimeline(FindChildTimeline(content.d[0], earthTimelineID), lifeTimelineID)";

            tl.left = ExecuteScriptGetNumber(script + ".left;");
            tl.right = ExecuteScriptGetNumber(script + ".right;");
            tl.height = ExecuteScriptGetNumber(script + ".height;");
            tl.realHeight = ExecuteScriptGetNumber(script + ".realHeight;");
            tl.y = ExecuteScriptGetNumber(script + ".y;");
            tl.realY = ExecuteScriptGetNumber(script + ".realY;");

            return tl;
        }

        public JsTimeline GetPrehistoryTimeline()
        {
            JsTimeline tl = new JsTimeline();
            string script = "return FindChildTimeline(FindChildTimeline(FindChildTimeline(content.d[0], earthTimelineID), " +
                            "lifeTimelineID), prehistoryTimelineID)";

            tl.left = ExecuteScriptGetNumber(script + ".left;");
            tl.right = ExecuteScriptGetNumber(script + ".right;");
            tl.height = ExecuteScriptGetNumber(script + ".height;");
            tl.realHeight = ExecuteScriptGetNumber(script + ".realHeight;");
            tl.y = ExecuteScriptGetNumber(script + ".y;");
            tl.realY = ExecuteScriptGetNumber(script + ".realY;");

            return tl;
        }

        public JsTimeline GetHumanityTimeline()
        {
            JsTimeline tl = new JsTimeline();
            string script = "return FindChildTimeline(FindChildTimeline(FindChildTimeline(FindChildTimeline(content.d[0], earthTimelineID), " +
                            "lifeTimelineID), prehistoryTimelineID), humanityTimelineID)";

            tl.left = ExecuteScriptGetNumber(script + ".left;");
            tl.right = ExecuteScriptGetNumber(script + ".right;");
            tl.height = ExecuteScriptGetNumber(script + ".height;");
            tl.realHeight = ExecuteScriptGetNumber(script + ".realHeight;");
            tl.y = ExecuteScriptGetNumber(script + ".y;");
            tl.realY = ExecuteScriptGetNumber(script + ".realY;");

            return tl;
        }

        // NOTE: Mark this test with TestMethod attribute, if need
        // to check error of mouse move actions in Selenium.
        public void TestRegimes_SeleniumMouseMoveError()
        {
            Point linkPosition = vcPageObj.CosmosLink.Location;

            for (int j = 5; j < 10; ++j)
            {
                for (int i = 35; i < 50; ++i)
                {
                    JsVisible v = vcPageObj.GetViewport();
                    action.MoveByOffset(linkPosition.X + i, linkPosition.Y + j).Click().Perform();
                    if (v != vcPageObj.GetViewport())
                    {
                        Console.WriteLine(i + " " + j);
                    }
                    action.SetDefault();
                }
            }
        }
    }

    [TestClass]
    [WebDriverSettings(BrowserType.Firefox)]
    public class RegimesTests_Firefox : RegimesTests
    {
        [ClassCleanup]
        public static void ClassCleanup()
        {
            Stop();
        }
    }

    [TestClass]
    [WebDriverSettings(BrowserType.InternetExplorer)]
    public class RegimesTests_IE : RegimesTests
    {
        [ClassCleanup]
        public static void ClassCleanup()
        {
            Stop();
        }
    }
}
